mod cli;

use std::process::ExitCode;

use clap::Parser;
use rsomics_common::Tool;

use cli::Cli;

fn main() -> ExitCode {
    Cli::parse().run()
}
