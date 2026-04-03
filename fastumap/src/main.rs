use anyhow::Result;
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "fastumap", about = "GPU-accelerated UMAP with cosine similarity")]
struct Cli {
    #[arg(short, long)]
    input: PathBuf,

    #[arg(short, long)]
    output: PathBuf,

    #[arg(long, default_value = "15")]
    n_neighbors: usize,

    #[arg(long, default_value = "0.1")]
    min_dist: f32,

    #[arg(long, default_value = "500")]
    n_epochs: usize,

    #[arg(long, default_value = "1.0")]
    learning_rate: f32,

    #[arg(long, default_value = "42")]
    seed: u64,

    #[arg(long, default_value = "0")]
    gpu_device: usize,

    #[arg(long, default_value = "false")]
    random_init: bool,

    /// Feature transform: none, sign, log
    #[arg(long, default_value = "none")]
    transform: String,

    /// Distance metric: cosine, hamming
    #[arg(long, default_value = "cosine")]
    metric: String,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    log::info!("Reading input: {}", cli.input.display());
    let (names, data) = fastumap::csv_io::read_csv(&cli.input, &cli.transform)?;

    let params = fastumap::UmapParams {
        n_neighbors: cli.n_neighbors,
        min_dist: cli.min_dist,
        n_epochs: cli.n_epochs,
        learning_rate: cli.learning_rate,
        seed: cli.seed,
        gpu_device: cli.gpu_device,
        random_init: cli.random_init,
        metric: cli.metric,
        transform: cli.transform,
        ..Default::default()
    };

    let result = fastumap::run_umap(&data, &params)?;

    log::info!("Writing output: {}", cli.output.display());
    fastumap::csv_io::write_csv(&cli.output, &names, &result.embedding)?;

    log::info!("Done!");
    Ok(())
}
