use anyhow::Result;
use clap::Args;

#[derive(Args)]
pub struct RsiArgs {
    pub input: String,
    #[arg(short = 'e', long = "explicit")]
    pub explicit: bool,
    #[arg(short = 'o', long = "output")]
    pub output: Option<String>,
}

pub fn run(_args: RsiArgs) -> Result<()> {
    todo!()
}
