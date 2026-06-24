use anyhow::Result;
use clap::Args;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(Args)]
pub struct ConvertArgs {
    #[arg(short = 'i', long = "in-fasta")]
    pub in_fasta: String,
    #[arg(long = "min-length", default_value_t = 0)]
    pub min_length: usize,
    #[arg(long = "max-length")]
    pub max_length: Option<usize>,
    #[arg(short = 'u', long = "usearch")]
    pub usearch: bool,
    #[arg(short = 's', long = "sample-fastas")]
    pub sample_fastas: bool,
}

pub fn run(args: ConvertArgs) -> Result<()> {
    let reader = BufReader::new(File::open(&args.in_fasta)?);
    let out_name = if args.usearch {
        "FilteredReads.forusearch.fna"
    } else {
        "FilteredReads.forsumaclust.fna"
    };
    let writer = BufWriter::new(File::create(out_name)?);

    let sample_dir = if args.sample_fastas {
        fs::create_dir_all("SampleFastas")?;
        Some("SampleFastas")
    } else {
        None
    };

    process(reader, writer, args.min_length, args.max_length, args.usearch, sample_dir)
}

fn process<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    min_length: usize,
    max_length: Option<usize>,
    usearch: bool,
    sample_dir: Option<&str>,
) -> Result<()> {
    let mut sample_handles: HashMap<String, BufWriter<File>> = HashMap::new();
    let mut counter: u64 = 1;

    let mut lines = reader.lines();
    while let Some(hdr_line) = lines.next() {
        let hdr_line = hdr_line?;
        if !hdr_line.starts_with('>') {
            continue;
        }
        let seq_line = match lines.next() {
            Some(l) => l?,
            None => break,
        };
        let seq = seq_line.trim_end();

        // Parse: >Sample TagPair counts_underscore
        let toks: Vec<&str> = hdr_line.splitn(4, ' ').collect();
        if toks.len() < 3 {
            continue;
        }
        let sample = &toks[0][1..];
        let size: u64 = toks[2]
            .split('_')
            .filter_map(|x| x.parse::<u64>().ok())
            .sum();

        if seq.len() < min_length {
            continue;
        }
        if let Some(max_len) = max_length {
            if seq.len() > max_len {
                continue;
            }
        }

        let out_hdr: String;
        let out_seq: String;
        if usearch {
            out_hdr = format!(">{};size={}", sample, size);
            out_seq = if let Some(max_len) = max_length {
                let mut s = seq.to_string();
                while s.len() < max_len {
                    s.push('N');
                }
                s
            } else {
                seq.to_string()
            };
        } else {
            out_hdr = format!(">{}:{} count={}", sample, counter, size);
            out_seq = seq.to_string();
            counter += 1;
        }

        writeln!(writer, "{}", out_hdr)?;
        writeln!(writer, "{}", out_seq)?;

        if let Some(dir) = sample_dir {
            if !sample_handles.contains_key(sample) {
                let path = format!("{}/{}.fixed.fasta", dir, sample);
                let fh = BufWriter::new(File::create(path)?);
                sample_handles.insert(sample.to_string(), fh);
            }
            let fh = sample_handles.get_mut(sample).unwrap();
            writeln!(fh, "{}", out_hdr)?;
            writeln!(fh, "{}", out_seq)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SMALL_FNA: &str = "\
>Sample1 Tag1-Tag2.Tag3-Tag4_1 5_4\n\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n\
>Sample1 Tag1-Tag2.Tag3-Tag4_2 3_2\n\
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n\
>Sample2 Tag5-Tag6.Tag7-Tag8_1 10_8\n\
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n\
>Sample2 Tag5-Tag6.Tag7-Tag8_2 1_0\n\
AAAA\n";

    #[test]
    fn test_sumaclust_output() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, None, false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], ">Sample1:1 count=9");
        assert_eq!(lines[1], "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
        assert_eq!(lines[2], ">Sample1:2 count=5");
        assert_eq!(lines[4], ">Sample2:3 count=18");
        assert_eq!(lines[6], ">Sample2:4 count=1");
        assert_eq!(lines[7], "AAAA");
        assert_eq!(lines.len(), 8);
    }

    #[test]
    fn test_usearch_output() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, None, true, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], ">Sample1;size=9");
        assert_eq!(lines[1], "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
        assert_eq!(lines[6], ">Sample2;size=1");
        assert_eq!(lines[7], "AAAA");
        assert_eq!(lines.len(), 8);
    }

    #[test]
    fn test_usearch_padding() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(65), true, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[1].len(), 65);
        assert!(lines[1].ends_with("NNNNN"));
        assert_eq!(lines[7].len(), 65);
        assert!(lines[7].starts_with("AAAA"));
        assert!(lines[7].ends_with('N'));
    }

    #[test]
    fn test_no_padding_in_sumaclust_mode() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(65), false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[1].len(), 60);
    }

    #[test]
    fn test_min_length_filter() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 10, None, false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 6);
        assert!(!s.contains("Sample2:4"));
        assert!(!s.contains("AAAA"));
    }

    #[test]
    fn test_max_length_filter() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(10), false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], ">Sample2:1 count=1");
        assert_eq!(lines[1], "AAAA");
    }
}
