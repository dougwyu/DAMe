use anyhow::Result;
use clap::Args;

#[derive(Args)]
pub struct SortArgs {
    #[arg(long = "fq")]
    pub fastq: String,
    #[arg(long = "primers")]
    pub primers: String,
    #[arg(long = "tags")]
    pub tags: String,
    #[arg(long = "keep-primers-seq")]
    pub keep_primers_seq: bool,
}

pub fn run(_args: SortArgs) -> Result<()> {
    todo!()
}
