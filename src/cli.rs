use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};

use rsomics_sc_pearson_residuals::{PearsonParams, open_output, run};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-sc-pearson-residuals",
    version,
    about = "Analytic Pearson residuals normalization of a single-cell count matrix (Lause 2021).",
    long_about = None
)]
pub struct Cli {
    /// 10x MTX directory (matrix.mtx[.gz], genes×cells).
    pub input: PathBuf,

    /// Output path for the dense residual matrix (MM array format); '-' for stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    /// Negative-binomial overdispersion parameter (shared across genes).
    #[arg(long, default_value_t = 100.0)]
    theta: f64,

    /// Symmetric clip bound. Defaults to sqrt(n_cells) per scanpy convention.
    #[arg(long)]
    clip: Option<f64>,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        if self.common.json {
            return Err(RsomicsError::InvalidInput(
                "--json is not supported: the residual matrix is written in MatrixMarket array \
                 format to --output (default stdout), the round-trippable single-cell convention"
                    .into(),
            ));
        }
        if self.theta <= 0.0 {
            return Err(RsomicsError::InvalidInput("--theta must be > 0".into()));
        }
        if let Some(c) = self.clip
            && c < 0.0
        {
            return Err(RsomicsError::InvalidInput("--clip must be >= 0".into()));
        }
        let params = PearsonParams {
            theta: self.theta,
            clip: self.clip,
        };
        let out = open_output(&self.output)?;
        let (genes, cells) = run(&self.input, &params, out)?;
        if !self.common.quiet {
            eprintln!(
                "pearson residuals: {cells} cells × {genes} genes (theta={}, clip={})",
                self.theta,
                self.clip.map_or_else(
                    || format!("sqrt({cells})={:.4}", (cells as f64).sqrt()),
                    |c| format!("{c}")
                )
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
