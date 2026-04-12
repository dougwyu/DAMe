use anyhow::Result;
use clap::Args;

#[derive(Args)]
pub struct DecollapseArgs {
    #[arg(long = "input")]
    pub input: String,
    #[arg(long = "out-fas", default_value = "Decollapsed.fasta")]
    pub out_fas: String,
}

pub fn run(_args: DecollapseArgs) -> Result<()> {
    todo!()
}
