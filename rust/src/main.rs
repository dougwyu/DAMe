use anyhow::Result;
use clap::{Parser, Subcommand};
use dame::{chimera_check, decollapse, filter, rsi, sort};

#[derive(Parser)]
#[command(name = "dame", about = "DNA Metabarcoding toolkit")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Sort(sort::SortArgs),
    Chimera(chimera_check::ChimeraArgs),
    Filter(filter::FilterArgs),
    Decollapse(decollapse::DecollapseArgs),
    Rsi(rsi::RsiArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Sort(args) => sort::run(args),
        Commands::Chimera(args) => chimera_check::run(args),
        Commands::Filter(args) => filter::run(args),
        Commands::Decollapse(args) => decollapse::run(args),
        Commands::Rsi(args) => rsi::run(args),
    }
}
