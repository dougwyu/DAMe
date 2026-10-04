//! `convert --per-pcr`: one FASTA record per PCR replicate, plus PCRinfo.txt.
//!
//! Input is Comparisons_<X>PCRs.fasta (or FilteredReads.fna) from `dame filter`.
//! Output reads keep their PCR as `sample=<sample>_PCR<k>` so a mapper such as
//! `vsearch --usearch_global --otutabout` builds a per-PCR table. The Python
//! implementation (python/dame/perpcr.py) must produce byte-identical output
//! and identical messages.

use anyhow::{bail, Context, Result};
use indexmap::IndexMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::filter::ps_info_rows;

pub const FASTA_OUT: &str = "FilteredReads.perpcr.fna";
pub const INFO_OUT: &str = "PCRinfo.txt";

pub const FILTERED_READS_WARNING: &str = "Warning: --per-pcr input looks like FilteredReads output, which has already been\nfiltered by --y/--t, so sequences that failed in a sample will appear as zeros.\nFor per-PCR OTU tables, use Comparisons_<X>PCRs.fasta instead.";

pub const USEARCH_NOTE: &str = "Note: -u has no effect with --per-pcr; per-PCR output always uses the\n;size=N;sample=<pcr_id>; label and is never padded.";

/// Returns the message for an invalid flag combination, or None.
pub fn flag_error(per_pcr: bool, ps_info: Option<&str>, sample_fastas: bool) -> Option<&'static str> {
    if per_pcr && sample_fastas {
        return Some("--sample-fastas cannot be used with --per-pcr");
    }
    if ps_info.is_some() && !per_pcr {
        return Some("--ps-info requires --per-pcr");
    }
    None
}

struct SampleState {
    tags: Vec<String>,
    reads: Vec<u64>,
}

/// PSinfo rows per sample in first-seen order; element k-1 is PCR k's (tag pair, pool).
type PsRows = IndexMap<String, Vec<(String, String)>>;

fn load_ps_info(ps_info: &str, x: usize) -> Result<PsRows> {
    let mut by_sample: IndexMap<String, IndexMap<usize, (String, String)>> = IndexMap::new();
    for (sample, pcr, tag_pair, pool) in ps_info_rows(ps_info, x)? {
        let rows = by_sample.entry(sample.clone()).or_default();
        if rows.contains_key(&pcr) {
            bail!("sample {} does not have exactly one PSinfo row for each PCR 1..{}", sample, x);
        }
        rows.insert(pcr, (tag_pair, pool));
    }
    let mut out = PsRows::new();
    for (sample, mut rows) in by_sample {
        if rows.len() != x || !(1..=x).all(|k| rows.contains_key(&k)) {
            bail!("sample {} does not have exactly one PSinfo row for each PCR 1..{}", sample, x);
        }
        let ordered = (1..=x).map(|k| rows.swap_remove(&k).unwrap()).collect();
        out.insert(sample, ordered);
    }
    Ok(out)
}

fn write_outputs(
    in_fasta: &str,
    ps_info: Option<&str>,
    min_length: usize,
    max_length: Option<usize>,
    fasta_tmp: &Path,
    info_tmp: &Path,
) -> Result<()> {
    let reader =
        BufReader::new(File::open(in_fasta).with_context(|| format!("opening {}", in_fasta))?);
    let mut out = BufWriter::new(
        File::create(fasta_tmp).with_context(|| format!("creating {}", fasta_tmp.display()))?,
    );

    let mut x: Option<usize> = None;
    let mut ps_rows: Option<PsRows> = None;
    let mut samples: IndexMap<String, SampleState> = IndexMap::new();
    let mut counter: u64 = 0;
    let mut ordinal: usize = 0;
    // (ordinal, sample, tag pairs, count parts) of a header awaiting its sequence line
    let mut pending: Option<(usize, String, Vec<String>, Vec<String>)> = None;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            ordinal += 1;
            pending = None;
            let toks: Vec<&str> = line.split_whitespace().collect();
            if toks.len() >= 3 {
                let tag_field = toks[1].rsplit_once('_').map_or(toks[1], |(head, _)| head);
                let tags = tag_field.split('.').map(String::from).collect();
                let parts = toks[2].split('_').map(String::from).collect();
                pending = Some((ordinal, toks[0][1..].to_string(), tags, parts));
            }
            continue;
        }
        let Some((ord, sample, tags, parts)) = pending.take() else {
            continue;
        };
        let seq = line.trim_end();

        let mut counts: Vec<i64> = Vec::with_capacity(parts.len());
        for part in &parts {
            match part.parse::<i64>() {
                Ok(c) => counts.push(c),
                Err(_) => bail!(
                    "record {} for sample {} has a non-integer count '{}'",
                    ord, sample, part
                ),
            }
        }
        if counts.len() != tags.len() {
            bail!(
                "record {} for sample {} has {} counts but {} tag pairs",
                ord, sample, counts.len(), tags.len()
            );
        }
        let xv = match x {
            None => {
                x = Some(tags.len());
                if let Some(p) = ps_info {
                    ps_rows = Some(load_ps_info(p, tags.len())?);
                }
                tags.len()
            }
            Some(xv) => {
                if tags.len() != xv {
                    bail!(
                        "record {} for sample {} has {} PCRs but earlier records have {}",
                        ord, sample, tags.len(), xv
                    );
                }
                xv
            }
        };

        if !samples.contains_key(&sample) {
            if let Some(rows) = &ps_rows {
                let Some(ps) = rows.get(&sample) else {
                    bail!(
                        "sample {} is in the input but not in {}",
                        sample,
                        ps_info.unwrap_or("")
                    );
                };
                for (k, (tag_pair, _pool)) in ps.iter().enumerate() {
                    if tags[k] != "empty-empty" && &tags[k] != tag_pair {
                        bail!(
                            "sample {} PCR {}: tag pair {} in the input but {} in PSinfo",
                            sample, k + 1, tags[k], tag_pair
                        );
                    }
                }
            }
            samples.insert(sample.clone(), SampleState { tags: tags.clone(), reads: vec![0; xv] });
        } else if samples[&sample].tags != tags {
            bail!(
                "sample {} has tag pairs {} in record {} but {} in an earlier record",
                sample, tags.join("."), ord, samples[&sample].tags.join(".")
            );
        }

        if seq.len() < min_length || max_length.map_or(false, |m| seq.len() > m) {
            continue;
        }
        let state = samples.get_mut(&sample).unwrap();
        for (k, &c) in counts.iter().enumerate() {
            if c > 0 {
                counter += 1;
                writeln!(
                    out,
                    ">{s}_PCR{k}.{n};size={c};sample={s}_PCR{k};",
                    s = sample, k = k + 1, n = counter, c = c
                )?;
                writeln!(out, "{}", seq)?;
                state.reads[k] += c as u64;
            }
        }
    }
    out.flush()?;

    let Some(x) = x else {
        bail!("no records in {}", in_fasta);
    };

    let mut info = BufWriter::new(
        File::create(info_tmp).with_context(|| format!("creating {}", info_tmp.display()))?,
    );
    match &ps_rows {
        None => {
            writeln!(info, "pcr_id\tsample\tpcr\ttag_pair\treads_pre_mapping")?;
            for (sample, state) in &samples {
                for k in 0..x {
                    let tag_pair = if state.tags[k] == "empty-empty" { "empty" } else { state.tags[k].as_str() };
                    writeln!(info, "{s}_PCR{k}\t{s}\t{k}\t{t}\t{r}", s = sample, k = k + 1, t = tag_pair, r = state.reads[k])?;
                }
            }
        }
        Some(rows) => {
            writeln!(info, "pcr_id\tsample\tpcr\ttag_pair\tpool\treads_pre_mapping")?;
            for (sample, ps) in rows {
                for (k, (tag_pair, pool)) in ps.iter().enumerate() {
                    let r = samples.get(sample).map_or(0, |s| s.reads[k]);
                    writeln!(info, "{s}_PCR{k}\t{s}\t{k}\t{t}\t{p}\t{r}", s = sample, k = k + 1, t = tag_pair, p = pool, r = r)?;
                }
            }
        }
    }
    info.flush()?;
    Ok(())
}

/// Writes FilteredReads.perpcr.fna and PCRinfo.txt into `out_dir`.
///
/// Outputs are written under `.tmp` names and renamed only when the whole
/// input has been read without error.
pub fn run_per_pcr(
    in_fasta: &str,
    ps_info: Option<&str>,
    min_length: usize,
    max_length: Option<usize>,
    usearch: bool,
    out_dir: &Path,
) -> Result<()> {
    let looks_filtered = Path::new(in_fasta)
        .file_name()
        .map_or(false, |n| n.to_string_lossy().starts_with("FilteredReads"));
    if looks_filtered {
        eprintln!("{}", FILTERED_READS_WARNING);
    }
    if usearch {
        eprintln!("{}", USEARCH_NOTE);
    }

    let fasta_out = out_dir.join(FASTA_OUT);
    let info_out = out_dir.join(INFO_OUT);
    let fasta_tmp = out_dir.join(format!("{}.tmp", FASTA_OUT));
    let info_tmp = out_dir.join(format!("{}.tmp", INFO_OUT));
    if let Err(e) = write_outputs(in_fasta, ps_info, min_length, max_length, &fasta_tmp, &info_tmp) {
        let _ = fs::remove_file(&fasta_tmp);
        let _ = fs::remove_file(&info_tmp);
        return Err(e);
    }
    fs::rename(&fasta_tmp, &fasta_out)?;
    fs::rename(&info_tmp, &info_out)?;
    Ok(())
}
