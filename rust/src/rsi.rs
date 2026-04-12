use anyhow::{Context, Result};
use clap::Args;
use ndarray::{Array2, Axis};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[derive(Args)]
pub struct RsiArgs {
    pub input: String,
    #[arg(short = 'e', long = "explicit")]
    pub explicit: bool,
    #[arg(short = 'o', long = "output")]
    pub output: Option<String>,
}

/// Compute the Renkonen Similarity Index between two replicates represented
/// as columns in `matrix`.
///
/// `matrix` is an N×2 array of integer counts (rows = sequences, cols = replicates).
/// `j` is the sample name (used only for diagnostic messages).
/// `a` and `b` are the replicate labels (used only for diagnostic messages).
///
/// Returns a value in [0, 1] where 0 = identical composition and 1 = no overlap.
pub fn compare(matrix: &Array2<i64>, j: &str, a: usize, b: usize) -> f64 {
    println!("Comparing replicates {a:?} and {b:?} from sample {j}");

    // Column sums (sum along rows → result has one value per column)
    let totals = matrix.sum_axis(Axis(0));
    let mut total_a = totals[0];
    let mut total_b = totals[1];

    if total_a == 0 {
        println!("This sample gave zero in replicate: {a:?}");
        total_a = 1;
    }
    if total_b == 0 {
        println!("This sample gave zero in replicate: {b:?}");
        total_b = 1;
    }

    let min_sum: f64 = matrix
        .rows()
        .into_iter()
        .map(|row| {
            let pa = row[0] as f64 / total_a as f64;
            let pb = row[1] as f64 / total_b as f64;
            pa.min(pb)
        })
        .sum();

    1.0 - min_sum
}

pub fn run(args: RsiArgs) -> Result<()> {
    // Read input file into rows of string fields
    let reader = BufReader::new(
        File::open(&args.input).with_context(|| format!("opening {}", args.input))?,
    );
    let data: Vec<Vec<String>> = reader
        .lines()
        .map(|l| {
            l.map(|s| s.split('\t').map(|f| f.to_string()).collect::<Vec<_>>())
        })
        .collect::<std::io::Result<Vec<_>>>()?;

    if data.is_empty() {
        println!("Input file is empty.");
        return Ok(());
    }

    let num_cols = data[0].len();
    // Format: sample | F0-R0 | count0 | F1-R1 | count1 | ... | seq
    // no_rep = (num_cols - 2) / 2
    let no_rep = (num_cols - 2) / 2;
    if no_rep < 2 {
        println!("There are no replicates in the file.");
        return Ok(());
    }

    // Group rows by sample name
    let mut sample_map: HashMap<String, Vec<Vec<String>>> = HashMap::new();
    for row in &data {
        sample_map
            .entry(row[0].clone())
            .or_default()
            .push(row.clone());
    }

    // Sort sample names for deterministic output
    let mut names: Vec<String> = sample_map.keys().cloned().collect();
    names.sort();

    // Result rows: for explicit mode (sample, rep_a, rep_b, rsi), else (sample, avg_rsi)
    let mut rkn_explicit: Vec<(String, usize, usize, f64)> = Vec::new();
    let mut rkn_avg: Vec<(String, f64)> = Vec::new();

    for name in &names {
        let subset = &sample_map[name];

        if args.explicit {
            for rep_a in 1..no_rep {
                for rep_b in (rep_a + 1)..=no_rep {
                    let matrix = build_matrix(subset, rep_a, rep_b);
                    let rsi = compare(&matrix, name, rep_a, rep_b);
                    rkn_explicit.push((name.clone(), rep_a, rep_b, rsi));
                }
            }
        } else {
            let mut total_rsi = 0.0_f64;
            let mut count = 0usize;
            for rep_a in 1..no_rep {
                for rep_b in (rep_a + 1)..=no_rep {
                    let matrix = build_matrix(subset, rep_a, rep_b);
                    total_rsi += compare(&matrix, name, rep_a, rep_b);
                    count += 1;
                }
            }
            rkn_avg.push((name.clone(), total_rsi / count as f64));
        }
    }

    let outfile = args
        .output
        .as_deref()
        .unwrap_or("RSI_output.txt");
    let mut out = File::create(outfile).with_context(|| format!("creating {outfile}"))?;

    if args.explicit {
        writeln!(out, "Sample\tReplicateA\tReplicateB\tRSI")?;
        for (sample, rep_a, rep_b, rsi) in &rkn_explicit {
            writeln!(out, "{sample}\t{rep_a}\t{rep_b}\t{rsi}")?;
        }
    } else {
        writeln!(out, "Sample\tRSI")?;
        for (sample, rsi) in &rkn_avg {
            writeln!(out, "{sample}\t{rsi}")?;
        }
    }

    Ok(())
}

/// Build an N×2 matrix from `subset` rows using columns `rep_a*2` and `rep_b*2`
/// (1-indexed replicate numbering, so col index = rep * 2).
fn build_matrix(subset: &[Vec<String>], rep_a: usize, rep_b: usize) -> Array2<i64> {
    let col_a = rep_a * 2;
    let col_b = rep_b * 2;
    let nrows = subset.len();
    let mut matrix = Array2::<i64>::zeros((nrows, 2));
    for (i, row) in subset.iter().enumerate() {
        matrix[[i, 0]] = row.get(col_a).and_then(|s| s.parse().ok()).unwrap_or(0);
        matrix[[i, 1]] = row.get(col_b).and_then(|s| s.parse().ok()).unwrap_or(0);
    }
    matrix
}
