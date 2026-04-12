use anyhow::Result;
use clap::Args;

#[derive(Args)]
pub struct ChimeraArgs {
    #[arg(long = "ps-info")]
    pub ps_info: String,
    #[arg(long = "x")]
    pub x: usize,
    #[arg(long = "p", default_value = "1")]
    pub p: usize,
}

pub fn run(_args: ChimeraArgs) -> Result<()> {
    todo!()
}
