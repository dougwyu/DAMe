use anyhow::{Context, Result};
use clap::Args;
use ahash::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[derive(Args)]
pub struct FilterArgs {
    /// Text file with tag combination info per PCR reaction per sample
    #[arg(long = "ps-info")]
    pub ps_info: String,
    /// Number of PCR rxns performed per sample
    #[arg(long = "x", default_value = "2")]
    pub x: usize,
    /// Number of PCR rxns sequence must be present in
    #[arg(long = "y", default_value = "1")]
    pub y: usize,
    /// Accepted for backward compatibility but ignored; the pool number is read from column 4 of PSinfo
    #[arg(long = "p", default_value = "1")]
    pub p: usize,
    /// Minimum count per unique sequence
    #[arg(long = "t", default_value = "1")]
    pub t: u32,
    /// Minimum sequence length
    #[arg(long = "l", default_value = "100")]
    pub l: usize,
    /// Use chimera-checked sorted collapsed files
    #[arg(long = "chimera-checked")]
    pub chimera_checked: bool,
}

/// Creates PS{n}_files.txt files (one per PCR replicate) from a PSinfo file.
/// When not chimera_checked: `pool{pool_num}/{tag1}_{tag2}.txt`
/// When chimera_checked: `{tag1}_{tag2}_{pool}.noChim.txt`
pub fn make_ps_num_files(ps_info: &str, x: usize, _p: usize, chimera_checked: bool) -> Result<()> {
    let mut outs: Vec<File> = (1..=x)
        .map(|i| {
            File::create(format!("PS{}_files.txt", i))
                .with_context(|| format!("creating PS{}_files.txt", i))
        })
        .collect::<Result<Vec<_>>>()?;

    let reader =
        BufReader::new(File::open(ps_info).with_context(|| format!("opening {}", ps_info))?);

    for (nr, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let nr = nr + 1; // 1-indexed
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let tag1 = parts[1];
        let tag2 = parts[2];
        let pool_num = parts[3];
        let residue = nr % x;
        let idx = if residue != 0 { residue - 1 } else { x - 1 };
        if !chimera_checked {
            writeln!(outs[idx], "pool{}/{}_{}.txt", pool_num, tag1, tag2)?;
        } else {
            writeln!(outs[idx], "{}_{}_{}.noChim.txt", tag1, tag2, pool_num)?;
        }
    }
    Ok(())
}

/// Reads PS{n}_files.txt files and returns a map of 0-indexed PCR replicate index -> list of file paths.
pub fn read_ps_num_files(x: usize) -> Result<HashMap<usize, Vec<String>>> {
    let mut ps_ins_lines: HashMap<usize, Vec<String>> = HashMap::default();
    for i in 0..x {
        let filename = format!("PS{}_files.txt", i + 1);
        let reader = BufReader::new(
            File::open(&filename).with_context(|| format!("opening {}", filename))?,
        );
        let lines: Vec<String> = reader
            .lines()
            .map(|l| l.map(|s| s.trim().to_string()))
            .collect::<std::io::Result<Vec<_>>>()?;
        ps_ins_lines.insert(i, lines);
    }
    Ok(ps_ins_lines)
}

/// Returns deduplicated list of sample names from PSinfo file, in order of first occurrence.
pub fn make_sample_name_array(ps_info: &str) -> Result<Vec<String>> {
    let reader =
        BufReader::new(File::open(ps_info).with_context(|| format!("opening {}", ps_info))?);
    let mut sample_names: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let name = line.split_whitespace().next().unwrap_or("").to_string();
        if !sample_names.contains(&name) {
            sample_names.push(name);
        }
    }
    Ok(sample_names)
}

/// Reads haplotype data for sample `i` across all PCR replicates.
/// Returns a map of replicate index (0-based) -> list of rows (each row is a Vec<String>).
pub fn read_haps_for_a_sample(
    x: usize,
    ps_ins_lines: &HashMap<usize, Vec<String>>,
    i: usize,
) -> Result<HashMap<usize, Vec<Vec<String>>>> {
    let mut haps: HashMap<usize, Vec<Vec<String>>> = HashMap::default();
    for j in 0..x {
        let entry = haps.entry(j).or_default();
        if let Some(lines) = ps_ins_lines.get(&j) {
            if i >= lines.len() {
                continue;
            }
            let path = &lines[i];
            if path == "empty" {
                continue;
            }
            if !std::path::Path::new(path).exists() {
                continue;
            }
            let reader = BufReader::new(
                File::open(path).with_context(|| format!("opening {}", path))?,
            );
            for line in reader.lines() {
                let line = line?;
                let row: Vec<String> = line.split_whitespace().map(|s| s.to_string()).collect();
                entry.push(row);
            }
        }
    }
    Ok(haps)
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplicateIndex {
    pub forward_tag: String,
    pub reverse_tag: String,
    pub counts_by_sequence: HashMap<String, i64>,
}

pub type SampleIndex = HashMap<usize, ReplicateIndex>;

/// Index each non-empty replicate once, preserving the first duplicate row.
pub fn index_haps(x: usize, haps: &HashMap<usize, Vec<Vec<String>>>) -> SampleIndex {
    let mut index = SampleIndex::default();
    for j in 0..x {
        let Some(rows) = haps.get(&j).filter(|rows| !rows.is_empty()) else {
            continue;
        };
        let (forward_tag, reverse_tag) = if rows[0].len() > 2 {
            (rows[0][1].clone(), rows[0][2].clone())
        } else {
            (String::new(), String::new())
        };
        let mut counts_by_sequence = HashMap::default();
        for row in rows.iter().filter(|row| row.len() >= 5) {
            let count = row[3].parse::<i64>().unwrap_or(0);
            counts_by_sequence.entry(row[4].clone()).or_insert(count);
        }
        index.insert(
            j,
            ReplicateIndex {
                forward_tag,
                reverse_tag,
                counts_by_sequence,
            },
        );
    }
    index
}

/// Return every sequence present in at least one indexed replicate.
pub fn all_sequences(index: &SampleIndex) -> HashSet<&str> {
    index
        .values()
        .flat_map(|replicate| replicate.counts_by_sequence.keys().map(String::as_str))
        .collect()
}

/// Writes comparison output files for all sequences of sample `i`.
#[allow(clippy::too_many_arguments)]
pub fn make_comparison_file(
    x: usize,
    seqs_all: &HashSet<&str>,
    index: &SampleIndex,
    out: &mut dyn Write,
    out_thresh: &mut dyn Write,
    out_yx: &mut dyn Write,
    out_fas: &mut dyn Write,
    out_thresh_fas: &mut dyn Write,
    out_yx_fas: &mut dyn Write,
    out_thresh_len_fas: &mut dyn Write,
    y: usize,
    t: u32,
    l: usize,
    sample_name: &[String],
    i: usize,
) -> Result<()> {
    let mut id_num: usize = 1;

    // Sort sequences for deterministic output
    let mut seqs_sorted: Vec<&str> = seqs_all.iter().copied().collect();
    seqs_sorted.sort();

    for seq in seqs_sorted {
        let sample = &sample_name[i];
        let mut line = format!("{}\t", sample);
        let mut line_fas_ids = format!(">{}\t", sample);
        let mut line_fas_counts = "\t".to_string();
        let mut y_count = 0usize;
        let mut t_count = 0usize;

        for j in 0..x {
            if let Some(replicate) = index.get(&j) {
                let count = if let Some(&count) = replicate.counts_by_sequence.get(seq) {
                    y_count += 1;
                    if (count as u32) < t {
                        t_count += 1;
                    }
                    count
                } else {
                    0
                };

                let fwd = replicate.forward_tag.as_str();
                let rev = replicate.reverse_tag.as_str();

                line.push_str(&format!("{}-{}\t{}\t", fwd, rev, count));

                if j < x - 1 {
                    line_fas_ids.push_str(&format!("{}-{}.", fwd, rev));
                    line_fas_counts.push_str(&format!("{}_", count));
                } else {
                    line_fas_ids.push_str(&format!("{}-{}_{}\t", fwd, rev, id_num));
                    line_fas_counts.push_str(&format!("{}\n{}", count, seq));
                }
            } else {
                // empty replicate
                line.push_str("empty\t0\t");
                if j < x - 1 {
                    line_fas_ids.push_str("empty-empty.");
                    line_fas_counts.push_str("0_");
                } else {
                    line_fas_ids.push_str(&format!("empty-empty_{}\t", id_num));
                    line_fas_counts.push_str(&format!("0\n{}", seq));
                }
            }
        }

        line.push_str(&format!("{}\n", seq));
        let line_fas = format!("{}{}\n", line_fas_ids, line_fas_counts);

        out.write_all(line.as_bytes())?;
        out_fas.write_all(line_fas.as_bytes())?;

        if y_count >= y {
            out_yx.write_all(line.as_bytes())?;
            out_yx_fas.write_all(line_fas.as_bytes())?;
        }

        if (y_count as i64 - t_count as i64) >= y as i64 {
            out_thresh.write_all(line.as_bytes())?;
            out_thresh_fas.write_all(line_fas.as_bytes())?;
            if seq.len() >= l {
                out_thresh_len_fas.write_all(line_fas.as_bytes())?;
            }
        }

        id_num += 1;
    }

    Ok(())
}

pub fn run(args: FilterArgs) -> Result<()> {
    let ps_info = &args.ps_info;
    let x = args.x;
    let y = args.y;
    let p = args.p;
    let t = args.t;
    let l = args.l;
    let chimera_checked = args.chimera_checked;

    let mut out = File::create(format!("Comparisons_{}PCRs.txt", x))
        .with_context(|| format!("creating Comparisons_{}PCRs.txt", x))?;
    let mut out_yx = File::create(format!("Comparisons_{}outOf{}PCRs.txt", y, x))
        .with_context(|| format!("creating Comparisons_{}outOf{}PCRs.txt", y, x))?;
    let mut out_thresh = File::create(format!(
        "Comparisons_{}outOf{}PCRs.countsThreshold{}.txt",
        y, x, t
    ))
    .with_context(|| "creating threshold file".to_string())?;
    let mut out_fas = File::create(format!("Comparisons_{}PCRs.fasta", x))
        .with_context(|| format!("creating Comparisons_{}PCRs.fasta", x))?;
    let mut out_yx_fas = File::create(format!("FilteredReads_atLeast{}.fasta", y))
        .with_context(|| format!("creating FilteredReads_atLeast{}.fasta", y))?;
    let mut out_thresh_fas = File::create(format!("FilteredReads_atLeast{}.threshold.fasta", y))
        .with_context(|| format!("creating FilteredReads_atLeast{}.threshold.fasta", y))?;
    let mut out_thresh_len_fas = File::create("FilteredReads.fna")
        .with_context(|| "creating FilteredReads.fna")?;

    make_ps_num_files(ps_info, x, p, chimera_checked)?;
    let ps_ins_lines = read_ps_num_files(x)?;
    let sample_name = make_sample_name_array(ps_info)?;

    let num_samples = ps_ins_lines.get(&0).map(|v| v.len()).unwrap_or(0);
    for i in 0..num_samples {
        let haps = read_haps_for_a_sample(x, &ps_ins_lines, i)?;
        let index = index_haps(x, &haps);
        let seqs_all = all_sequences(&index);
        make_comparison_file(
            x,
            &seqs_all,
            &index,
            &mut out,
            &mut out_thresh,
            &mut out_yx,
            &mut out_fas,
            &mut out_thresh_fas,
            &mut out_yx_fas,
            &mut out_thresh_len_fas,
            y,
            t,
            l,
            &sample_name,
            i,
        )?;
    }

    Ok(())
}
