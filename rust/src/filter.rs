use anyhow::{Context, Result};
use clap::Args;
use ahash::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Filter amplicons across PCR replicates by abundance and sequence length.
///
/// Reads the per-tag-pair collapsed files produced by `dame sort`, groups them
/// by sample and PCR replicate, and writes seven output files at different
/// filtering stringencies (all seqs / passed-y / passed-y-and-t / length-filtered).
#[derive(Args)]
pub struct FilterArgs {
    /// PSinfo file: tab-separated columns (SampleName, FwdTagName, RevTagName, PoolNumber),
    /// one PCR reaction per line.  Same file used by `dame sort --ps-info`.
    #[arg(long = "ps-info")]
    pub ps_info: String,

    /// Number of PCR replicates per sample [default: 2]
    #[arg(long = "x", default_value = "2")]
    pub x: usize,

    /// Minimum number of replicates a sequence must appear in to pass [default: 1]
    #[arg(long = "y", default_value = "1")]
    pub y: usize,

    /// Number of sequencing pools [default: 1]
    #[arg(long = "p", default_value = "1")]
    pub p: usize,

    /// Minimum read count per sequence per replicate; counts below this threshold
    /// increment the below-threshold counter used in conjunction with --y [default: 1]
    #[arg(long = "t", default_value = "1")]
    pub t: u32,

    /// Minimum sequence length in nucleotides (applied to FilteredReads.fna only) [default: 100]
    #[arg(long = "l", default_value = "100")]
    pub l: usize,

    /// Use chimera-checked input files ({tag1}_{tag2}_{pool}.noChim.txt) instead of
    /// the default sort output (pool{n}/{tag1}_{tag2}.txt)
    #[arg(long = "chimera-checked")]
    pub chimera_checked: bool,

    /// Output directory; auto-named Filter_min{Y}PCRs_min{T}copies_{marker} if omitted,
    /// where marker is derived from the --ps-info filename
    /// (e.g. PCRsetsInfo_MIFISH.txt → Filter_min1PCRs_min1copies_MIFISH)
    #[arg(short = 'o', long = "outdir")]
    pub outdir: Option<String>,
}

/// Derives the auto output-directory name from filter parameters and the PSinfo filename.
/// Pattern: `Filter_min{y}PCRs_min{t}copies_{marker}`
/// where marker strips known extensions and the `PCRsetsInfo_` prefix.
pub fn auto_outdir(ps_info: &str, y: usize, t: u32) -> String {
    let basename = std::path::Path::new(ps_info)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(ps_info);
    let stem = basename
        .strip_suffix(".txt.gz")
        .or_else(|| basename.strip_suffix(".txt"))
        .or_else(|| basename.strip_suffix(".gz"))
        .unwrap_or(basename);
    let marker = stem.strip_prefix("PCRsetsInfo_").unwrap_or(stem);
    format!("Filter_min{}PCRs_min{}copies_{}", y, t, marker)
}

/// Creates PS{n}_files.txt files (one per PCR replicate) from a PSinfo file.
/// When not chimera_checked: `pool{pool_num}/{tag1}_{tag2}.txt`
/// When chimera_checked: `{tag1}_{tag2}_{pool}.noChim.txt`
pub fn make_ps_num_files(ps_info: &str, x: usize, _p: usize, chimera_checked: bool, outdir: &str) -> Result<()> {
    let mut outs: Vec<File> = (1..=x)
        .map(|i| {
            let path = format!("{}/PS{}_files.txt", outdir, i);
            File::create(&path).with_context(|| format!("creating {}", path))
        })
        .collect::<Result<Vec<_>>>()?;

    let reader =
        BufReader::new(File::open(ps_info).with_context(|| format!("opening {}", ps_info))?);

    // Group file entries by sample name, maintaining first-occurrence order
    let mut sample_order: Vec<String> = Vec::new();
    let mut sample_entries: std::collections::HashMap<String, Vec<String>> =
        std::collections::HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let sample = parts[0].to_string();
        let tag1 = parts[1];
        let tag2 = parts[2];
        let pool_num = parts[3];
        let entry = if !chimera_checked {
            format!("pool{}/{}_{}.txt", pool_num, tag1, tag2)
        } else {
            format!("{}_{}_{}.noChim.txt", tag1, tag2, pool_num)
        };
        if !sample_entries.contains_key(&sample) {
            sample_order.push(sample.clone());
            sample_entries.insert(sample.clone(), Vec::new());
        }
        sample_entries.get_mut(&sample).unwrap().push(entry);
    }

    // Write: k-th replicate of each sample → PSk+1
    for sample in &sample_order {
        let entries = &sample_entries[sample];
        for k in 0..x {
            if k < entries.len() {
                writeln!(outs[k], "{}", entries[k])?;
            } else {
                writeln!(outs[k], "empty")?;
            }
        }
    }
    Ok(())
}

/// Reads PS{n}_files.txt files and returns a map of 0-indexed PCR replicate index -> list of file paths.
pub fn read_ps_num_files(x: usize, outdir: &str) -> Result<HashMap<usize, Vec<String>>> {
    let mut ps_ins_lines: HashMap<usize, Vec<String>> = HashMap::default();
    for i in 0..x {
        let filename = format!("{}/PS{}_files.txt", outdir, i + 1);
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
                let row: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
                entry.push(row);
            }
        }
    }
    Ok(haps)
}

/// Return type for `get_seqs_sets_and_fr_counts`: (seqs_all, F, R, seq_counts)
pub type SeqsSetsAndFRCounts = (
    HashSet<String>,
    HashMap<usize, String>,
    HashMap<usize, String>,
    HashMap<usize, HashMap<String, String>>,
);

/// Builds the set of all unique sequences and per-replicate lookup maps from haplotype data.
///
/// Returns: (seqs_all, F, R, seq_counts)
/// - seqs_all:   HashSet of all unique sequences across all replicates
/// - F:          replicate index -> forward tag name (column 1 of first row)
/// - R:          replicate index -> reverse tag name (column 2 of first row)
/// - seq_counts: replicate index -> {sequence -> count_str}  (O(1) lookup at call site)
pub fn get_seqs_sets_and_fr_counts(
    x: usize,
    haps: &HashMap<usize, Vec<Vec<String>>>,
) -> SeqsSetsAndFRCounts {
    let mut f: HashMap<usize, String> = HashMap::default();
    let mut r: HashMap<usize, String> = HashMap::default();
    let mut seq_counts: HashMap<usize, HashMap<String, String>> = HashMap::default();
    let mut seqs_all: HashSet<String> = HashSet::default();

    for j in 0..x {
        if let Some(hap_j) = haps.get(&j) {
            if !hap_j.is_empty() {
                if hap_j[0].len() > 2 {
                    f.insert(j, hap_j[0][1].clone());
                    r.insert(j, hap_j[0][2].clone());
                }
                let mut j_map: HashMap<String, String> = HashMap::default();
                for row in hap_j {
                    if row.len() >= 5 {
                        seqs_all.insert(row[4].clone());
                        j_map.insert(row[4].clone(), row[3].clone());
                    }
                }
                seq_counts.insert(j, j_map);
            }
        }
    }

    (seqs_all, f, r, seq_counts)
}

/// Writes comparison output files for all sequences of sample `i`.
#[allow(clippy::too_many_arguments)]
pub fn make_comparison_file(
    x: usize,
    seqs_all: &HashSet<String>,
    f: &HashMap<usize, String>,
    r: &HashMap<usize, String>,
    seq_counts: &HashMap<usize, HashMap<String, String>>,
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
            let replicate_has_data = f.contains_key(&j);

            if replicate_has_data {
                // O(1) lookup: seq_counts[j] is a HashMap<seq -> count_str>
                let count_str = seq_counts
                    .get(&j)
                    .and_then(|m| m.get(*seq))
                    .map(|s| s.as_str())
                    .unwrap_or("0");

                let count: i64 = if count_str != "0" {
                    y_count += 1;
                    let c: i64 = count_str.parse().unwrap_or(0);
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

    let outdir = args.outdir.unwrap_or_else(|| auto_outdir(ps_info, y, t));
    std::fs::create_dir_all(&outdir)
        .with_context(|| format!("creating output directory {}", outdir))?;

    let out_path = |name: &str| format!("{}/{}", outdir, name);

    let mut out = File::create(out_path(&format!("Comparisons_{}PCRs.txt", x)))
        .with_context(|| format!("creating Comparisons_{}PCRs.txt", x))?;
    let mut out_yx = File::create(out_path(&format!("Comparisons_{}outOf{}PCRs.txt", y, x)))
        .with_context(|| format!("creating Comparisons_{}outOf{}PCRs.txt", y, x))?;
    let mut out_thresh = File::create(out_path(&format!(
        "Comparisons_{}outOf{}PCRs.countsThreshold{}.txt",
        y, x, t
    )))
    .with_context(|| "creating threshold file".to_string())?;
    let mut out_fas = File::create(out_path(&format!("Comparisons_{}PCRs.fasta", x)))
        .with_context(|| format!("creating Comparisons_{}PCRs.fasta", x))?;
    let mut out_yx_fas = File::create(out_path(&format!("FilteredReads_atLeast{}.fasta", y)))
        .with_context(|| format!("creating FilteredReads_atLeast{}.fasta", y))?;
    let mut out_thresh_fas =
        File::create(out_path(&format!("FilteredReads_atLeast{}.threshold.fasta", y)))
            .with_context(|| format!("creating FilteredReads_atLeast{}.threshold.fasta", y))?;
    let mut out_thresh_len_fas = File::create(out_path("FilteredReads.fna"))
        .with_context(|| "creating FilteredReads.fna")?;

    make_ps_num_files(ps_info, x, p, chimera_checked, &outdir)?;
    let ps_ins_lines = read_ps_num_files(x, &outdir)?;
    let sample_name = make_sample_name_array(ps_info)?;

    let num_samples = ps_ins_lines.get(&0).map(|v| v.len()).unwrap_or(0);
    for i in 0..num_samples {
        let haps = read_haps_for_a_sample(x, &ps_ins_lines, i)?;
        let (seqs_all, f, r, seq_counts) = get_seqs_sets_and_fr_counts(x, &haps);
        make_comparison_file(
            x,
            &seqs_all,
            &f,
            &r,
            &seq_counts,
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

    eprintln!("Output written to: {}", outdir);
    Ok(())
}
