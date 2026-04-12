use anyhow::Result;
use clap::Args;

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

pub fn run(_args: FilterArgs) -> Result<()> {
    todo!()
}
