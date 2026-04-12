use anyhow::Result;
use clap::Args;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(Args)]
pub struct DecollapseArgs {
    #[arg(long = "input")]
    pub input: String,
    #[arg(long = "out-fas", default_value = "Decollapsed.fasta")]
    pub out_fas: String,
}

pub fn run(args: DecollapseArgs) -> Result<()> {
    let in_file = BufReader::new(File::open(&args.input)?);
    let out_file = BufWriter::new(File::create(&args.out_fas)?);
    process(in_file, out_file)
}

fn process<R: BufRead, W: Write>(reader: R, mut writer: W) -> Result<()> {
    let mut seq_id: u64 = 0;
    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        // Columns: PrimerName, Tag1, Tag2, Frequency, Sequence
        let tag1 = fields[1];
        let tag2 = fields[2];
        let freq: u64 = fields[3].parse()?;
        let seq = fields[4];
        for _ in 0..freq {
            seq_id += 1;
            writeln!(
                writer,
                ">{tag1}.{tag2}.{freq}_{seq_id}\n{seq}",
            )?;
        }
    }
    Ok(())
}
