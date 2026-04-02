use anyhow::{Context, Result};
use ndarray::Array2;
use std::path::Path;

pub fn read_csv(path: &Path) -> Result<(Vec<String>, Array2<f32>)> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .with_context(|| format!("Failed to open {}", path.display()))?;

    let mut names: Vec<String> = Vec::new();
    let mut data: Vec<Vec<f32>> = Vec::new();

    for (i, result) in reader.records().enumerate() {
        let record = result.with_context(|| format!("Failed to read row {}", i))?;
        let name = record.get(0).unwrap_or("").to_string();
        names.push(name);

        // Skip column 1 (depth), take columns 2..end as features
        let features: Vec<f32> = record
            .iter()
            .skip(2)
            .map(|s| s.trim().parse::<f32>().unwrap_or(0.0))
            .collect();
        data.push(features);
    }

    let n_rows = data.len();
    let n_cols = data.first().map(|r| r.len()).unwrap_or(0);
    log::info!("Read {} cells with {} features", n_rows, n_cols);

    let flat: Vec<f32> = data.into_iter().flatten().collect();
    let matrix = Array2::from_shape_vec((n_rows, n_cols), flat)
        .context("Failed to create matrix from CSV data")?;

    Ok((names, matrix))
}

pub fn write_csv(path: &Path, names: &[String], coords: &[[f32; 2]]) -> Result<()> {
    let mut writer = csv::Writer::from_path(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;

    writer.write_record(["cell_name", "umap_x", "umap_y"])?;
    for (name, coord) in names.iter().zip(coords.iter()) {
        writer.write_record(&[name.as_str(), &coord[0].to_string(), &coord[1].to_string()])?;
    }
    writer.flush()?;
    Ok(())
}
