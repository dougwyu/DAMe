use anyhow::{anyhow, Context, Result};
use clap::Args;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::process::Command;

#[derive(Args)]
pub struct ChimeraArgs {
    #[arg(long = "ps-info")]
    pub ps_info: String,
    #[arg(long = "x")]
    pub x: usize,
    #[arg(long = "p", default_value = "1")]
    pub p: usize,
}

/// Creates PS{n}.tags.txt files (one per PCR replicate) from a PSinfo file.
/// Output format per line: `tag1\ttag2\n` (no pool number).
pub fn make_tag_files(ps_info: &str, x: usize) -> Result<()> {
    let mut outs: Vec<File> = (1..=x)
        .map(|i| {
            File::create(format!("PS{}.tags.txt", i))
                .with_context(|| format!("creating PS{}.tags.txt", i))
        })
        .collect::<Result<Vec<_>>>()?;

    let reader = BufReader::new(File::open(ps_info).with_context(|| format!("opening {}", ps_info))?);

    for (nr, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let nr = nr + 1; // 1-indexed
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let tag1 = parts[1];
        let tag2 = parts[2];
        let residue = nr % x;
        let idx = if residue != 0 { residue - 1 } else { x - 1 };
        writeln!(outs[idx], "{}\t{}", tag1, tag2)?;
    }
    Ok(())
}

/// Creates PS{n}.tags.txt files including pool number.
/// Output format per line: `tag1\ttag2\tpool\n`.
pub fn make_tag_files_with_pools(ps_info: &str, x: usize) -> Result<()> {
    let mut outs: Vec<File> = (1..=x)
        .map(|i| {
            File::create(format!("PS{}.tags.txt", i))
                .with_context(|| format!("creating PS{}.tags.txt", i))
        })
        .collect::<Result<Vec<_>>>()?;

    let reader = BufReader::new(File::open(ps_info).with_context(|| format!("opening {}", ps_info))?);

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
        let pool = parts[3];
        let residue = nr % x;
        let idx = if residue != 0 { residue - 1 } else { x - 1 };
        writeln!(outs[idx], "{}\t{}\t{}", tag1, tag2, pool)?;
    }
    Ok(())
}

/// Builds Pool{p}.fasta from the PS{n}.tags.txt files.
pub fn make_size_out_fastas(p: usize, x: usize) -> Result<()> {
    let mut outs: Vec<File> = (1..=p)
        .map(|pool| {
            File::create(format!("Pool{}.fasta", pool))
                .with_context(|| format!("creating Pool{}.fasta", pool))
        })
        .collect::<Result<Vec<_>>>()?;

    for num in 0..x {
        let tag_file = format!("PS{}.tags.txt", num + 1);
        let reader = BufReader::new(
            File::open(&tag_file).with_context(|| format!("opening {}", tag_file))?,
        );

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 2 {
                continue;
            }
            let tag1 = parts[0];
            let tag2 = parts[1];

            let hap_path = if p > 1 && parts.len() >= 3 {
                let pool_num = parts[2];
                format!("./pool{}/{}_{}.txt", pool_num, tag1, tag2)
            } else {
                format!("{}_{}.txt", tag1, tag2)
            };

            if !std::path::Path::new(&hap_path).exists() {
                continue;
            }

            let hap_reader = BufReader::new(
                File::open(&hap_path).with_context(|| format!("opening {}", hap_path))?,
            );

            let pool_idx = if p > 1 && parts.len() >= 3 {
                parts[2].parse::<usize>().unwrap_or(1) - 1
            } else {
                0
            };

            let mut id_num = 1usize;
            for seq_line in hap_reader.lines() {
                let seq_line = seq_line?;
                let seq_line = seq_line.trim();
                if seq_line.is_empty() {
                    continue;
                }
                let seq_parts: Vec<&str> = seq_line.split('\t').collect();
                if seq_parts.len() < 5 {
                    continue;
                }
                let primer = seq_parts[0];
                let freq = seq_parts[3];
                let sequence = seq_parts[4];
                let header = format!(
                    ">{}_{}_{}_{};size={}",
                    primer, tag1, tag2, id_num, freq
                );
                writeln!(outs[pool_idx], "{}", header)?;
                writeln!(outs[pool_idx], "{}", sequence)?;
                id_num += 1;
            }
        }
    }
    Ok(())
}

/// Runs usearch twice per pool (sortsize + uchime).
pub fn sort_fasta(p: usize) -> Result<()> {
    for pool in 1..=p {
        let input_file = format!("Pool{}.fasta", pool);
        let sort_output = format!("Pool{}.sort.fasta", pool);

        let status = Command::new("usearch")
            .args([
                "--sortsize",
                &input_file,
                "--output",
                &sort_output,
            ])
            .output()
            .map_err(|e| {
                if e.kind() == std::io::ErrorKind::NotFound {
                    anyhow!("usearch not found in PATH. Please install usearch and ensure it is accessible.")
                } else {
                    anyhow!("Failed to run usearch --sortsize: {}", e)
                }
            })?;

        fs::write(format!("sort{}.out", pool), &status.stdout)?;
        fs::write(format!("sort{}.err", pool), &status.stderr)?;

        let chim_output = format!("Pool{}.Chim.fasta", pool);
        let nochim_output = format!("Pool{}.noChim.fasta", pool);

        let status2 = Command::new("usearch")
            .args([
                "-uchime",
                &sort_output,
                "-chimeras",
                &chim_output,
                "-nonchimeras",
                &nochim_output,
            ])
            .output()
            .map_err(|e| {
                if e.kind() == std::io::ErrorKind::NotFound {
                    anyhow!("usearch not found in PATH. Please install usearch and ensure it is accessible.")
                } else {
                    anyhow!("Failed to run usearch -uchime: {}", e)
                }
            })?;

        fs::write(format!("chimeraCheck{}.out", pool), &status2.stdout)?;
        fs::write(format!("chimeraCheck{}.err", pool), &status2.stderr)?;
    }
    Ok(())
}

/// Reads Pool{p}.noChim.fasta and writes Pool{p}.noChim.oneLiner.fasta
/// where multi-line sequences are joined into a single line.
pub fn make_fas_seq_one_line(p: usize) -> Result<()> {
    for pool in 1..=p {
        let input_path = format!("Pool{}.noChim.fasta", pool);
        let output_path = format!("Pool{}.noChim.oneLiner.fasta", pool);

        let reader = BufReader::new(
            File::open(&input_path).with_context(|| format!("opening {}", input_path))?,
        );
        let mut out =
            File::create(&output_path).with_context(|| format!("creating {}", output_path))?;

        let mut seq = String::new();
        for line in reader.lines() {
            let line = line?;
            let line = line.trim_end_matches('\r');
            if line.starts_with('>') {
                if !seq.is_empty() {
                    writeln!(out, "{}", seq)?;
                    seq.clear();
                }
                writeln!(out, "{}", line)?;
            } else {
                seq.push_str(line);
            }
        }
        if !seq.is_empty() {
            writeln!(out, "{}", seq)?;
        }
    }
    Ok(())
}

type NoChimEntry = (Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>);

/// Parses Pool{p}.noChim.oneLiner.fasta and creates {tag1}_{tag2}_{pool}.noChim.txt files.
///
/// FASTA header format: `>{primer}_{tag1}_{tag2}_{idnum};size={freq}`
#[allow(unused_assignments)]
pub fn make_no_chim_haps(p: usize) -> Result<()> {
    // HAP: key = "tag1_tag2_pool" → (primers, tag1s, tag2s, freqs, seqs)
    let mut hap: HashMap<String, NoChimEntry> = HashMap::new();

    // Preserve insertion order for deterministic output
    let mut hap_order: Vec<String> = Vec::new();

    for pool in 1..=p {
        let input_path = format!("Pool{}.noChim.oneLiner.fasta", pool);
        let reader = BufReader::new(
            File::open(&input_path).with_context(|| format!("opening {}", input_path))?,
        );

        let mut primer_name = String::new();
        let mut tag_name1 = String::new();
        let mut tag_name2 = String::new();
        let mut freq = String::new();
        let mut tag_hap_key = String::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim_end_matches('\r');
            if let Some(rest) = line.strip_prefix('>') {
                // format: {primer}_{tag1}_{tag2}_{idnum};size={freq}
                let parts: Vec<&str> = rest.split('_').collect();
                if parts.len() >= 3 {
                    primer_name = parts[0].to_string();
                    tag_name1 = parts[1].to_string();
                    tag_name2 = parts[2].to_string();
                }
                // Extract freq from "...;size=N"
                freq = rest
                    .split('=')
                    .nth(1)
                    .unwrap_or("")
                    .to_string();
                tag_hap_key = format!("{}_{}_{}", tag_name1, tag_name2, pool);

                let entry = hap.entry(tag_hap_key.clone()).or_insert_with(|| {
                    hap_order.push(tag_hap_key.clone());
                    (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new())
                });
                entry.0.push(primer_name.clone());
                entry.1.push(tag_name1.clone());
                entry.2.push(tag_name2.clone());
                entry.3.push(freq.clone());
            } else if !tag_hap_key.is_empty() {
                if let Some(entry) = hap.get_mut(&tag_hap_key) {
                    entry.4.push(format!("{}\n", line));
                }
            }
        }
    }

    for key in &hap_order {
        let entry = &hap[key];
        let out_path = format!("{}.noChim.txt", key);
        let mut out = File::create(&out_path).with_context(|| format!("creating {}", out_path))?;
        for i in 0..entry.0.len() {
            if i < entry.4.len() {
                write!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    entry.0[i], entry.1[i], entry.2[i], entry.3[i], entry.4[i]
                )?;
            }
        }
    }

    Ok(())
}

pub fn run(args: ChimeraArgs) -> Result<()> {
    let ps_info = &args.ps_info;
    let x = args.x;
    let p = args.p;

    if p == 1 {
        make_tag_files(ps_info, x)?;
    } else {
        make_tag_files_with_pools(ps_info, x)?;
    }

    make_size_out_fastas(p, x)?;
    sort_fasta(p)?;
    make_fas_seq_one_line(p)?;
    make_no_chim_haps(p)?;

    Ok(())
}
