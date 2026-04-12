use anyhow::{Context, Result};
use clap::Args;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[derive(Args)]
pub struct FilterArgs {
    #[arg(long = "ps-info")]
    pub ps_info: String,
    #[arg(long = "x", default_value = "2")]
    pub x: usize,
    #[arg(long = "y", default_value = "1")]
    pub y: usize,
    #[arg(long = "p", default_value = "1")]
    pub p: usize,
    #[arg(long = "t", default_value = "1")]
    pub t: u32,
    #[arg(long = "l", default_value = "100")]
    pub l: usize,
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
        let parts: Vec<&str> = line.split('\t').collect();
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
    let mut ps_ins_lines: HashMap<usize, Vec<String>> = HashMap::new();
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
        let name = line.split('\t').next().unwrap_or("").to_string();
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
    let mut haps: HashMap<usize, Vec<Vec<String>>> = HashMap::new();
    for j in 0..x {
        haps.insert(j, Vec::new());
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
            let entry = haps.get_mut(&j).unwrap();
            for line in reader.lines() {
                let line = line?;
                let row: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
                entry.push(row);
            }
        }
    }
    Ok(haps)
}

/// Return type for `get_seqs_sets_and_fr_counts`: (seqs_all, F, R, counts, seqs)
pub type SeqsSetsAndFRCounts = (
    HashSet<String>,
    HashMap<usize, String>,
    HashMap<usize, String>,
    HashMap<usize, Vec<String>>,
    HashMap<usize, Vec<String>>,
);

/// Builds the set of all unique sequences and per-replicate F, R, counts, seqs from haplotype data.
///
/// Returns: (seqs_all, F, R, counts, seqs)
/// - seqs_all: HashSet of all unique sequences across all replicates
/// - F: map of replicate index -> forward tag (tag1) from first row
/// - R: map of replicate index -> reverse tag (tag2) from first row
/// - counts: map of replicate index -> list of count strings
/// - seqs: map of replicate index -> list of sequence strings
pub fn get_seqs_sets_and_fr_counts(
    x: usize,
    haps: &HashMap<usize, Vec<Vec<String>>>,
) -> SeqsSetsAndFRCounts {
    let mut f: HashMap<usize, String> = HashMap::new();
    let mut r: HashMap<usize, String> = HashMap::new();
    let mut counts: HashMap<usize, Vec<String>> = HashMap::new();
    let mut seqs: HashMap<usize, Vec<String>> = HashMap::new();
    let mut seqs_all: HashSet<String> = HashSet::new();

    for j in 0..x {
        if let Some(hap_j) = haps.get(&j) {
            if !hap_j.is_empty() {
                let mut j_seqs: Vec<String> = Vec::new();
                let mut j_counts: Vec<String> = Vec::new();
                // Get F and R from first row (columns 1 and 2)
                if hap_j[0].len() > 2 {
                    f.insert(j, hap_j[0][1].clone());
                    r.insert(j, hap_j[0][2].clone());
                }
                for row in hap_j {
                    if row.len() >= 5 {
                        j_counts.push(row[3].clone());
                        j_seqs.push(row[4].clone());
                        seqs_all.insert(row[4].clone());
                    }
                }
                counts.insert(j, j_counts);
                seqs.insert(j, j_seqs);
            }
        }
    }

    (seqs_all, f, r, counts, seqs)
}

/// Writes comparison output files for all sequences of sample `i`.
#[allow(clippy::too_many_arguments)]
pub fn make_comparison_file(
    x: usize,
    seqs_all: &HashSet<String>,
    haps: &HashMap<usize, Vec<Vec<String>>>,
    f: &HashMap<usize, String>,
    r: &HashMap<usize, String>,
    counts: &HashMap<usize, Vec<String>>,
    seqs: &HashMap<usize, Vec<String>>,
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
    let mut seqs_sorted: Vec<&String> = seqs_all.iter().collect();
    seqs_sorted.sort();

    for seq in &seqs_sorted {
        let sample = &sample_name[i];
        let mut line = format!("{}\t", sample);
        let mut line_fas_ids = format!(">{}\t", sample);
        let mut line_fas_counts = "\t".to_string();
        let mut y_count = 0usize;
        let mut t_count = 0usize;

        for j in 0..x {
            let hap_j_empty = haps.get(&j).map(|v| v.is_empty()).unwrap_or(true);

            if !hap_j_empty {
                // Find position of seq in this replicate's seqs list
                let pos = seqs
                    .get(&j)
                    .and_then(|sv| sv.iter().position(|s| s == *seq));

                let count: i64 = if let Some(p) = pos {
                    y_count += 1;
                    let c_str = counts
                        .get(&j)
                        .and_then(|cv| cv.get(p))
                        .map(|s| s.as_str())
                        .unwrap_or("0");
                    let c: i64 = c_str.parse().unwrap_or(0);
                    if (c as u32) < t {
                        t_count += 1;
                    }
                    c
                } else {
                    0
                };

                let fwd = f.get(&j).map(|s| s.as_str()).unwrap_or("");
                let rev = r.get(&j).map(|s| s.as_str()).unwrap_or("");

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
        let (seqs_all, f, r, counts, seqs) = get_seqs_sets_and_fr_counts(x, &haps);
        make_comparison_file(
            x,
            &seqs_all,
            &haps,
            &f,
            &r,
            &counts,
            &seqs,
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
