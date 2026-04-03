use anyhow::{Context, Result};
use cudarc::cublas::safe::GemmConfig;
use cudarc::cublas::sys;
use cudarc::cublas::Gemm;
use cudarc::driver::{LaunchConfig, PushKernelArg};
use rayon::prelude::*;

use super::context::GpuContext;
use crate::umap::knn_graph::KnnGraph;

pub fn compute_cosine_knn(
    gpu: &GpuContext,
    data: &[f32],
    n_rows: usize,
    n_cols: usize,
    k: usize,
    normalize: bool,
) -> Result<KnnGraph> {
    log::info!(
        "Computing cosine k-NN on GPU: {} x {}, k={}",
        n_rows,
        n_cols,
        k
    );

    // Upload data to GPU (row-major: n_rows x n_cols)
    let mut d_data = gpu
        .stream
        .clone_htod(data)
        .context("Failed to upload data to GPU")?;

    // Compute row norms (needed for output normalization)
    let mut d_norms = gpu.stream.alloc_zeros::<f32>(n_rows)?;
    {
        let block_size = 256u32;
        let grid_size = ((n_rows as u32) + block_size - 1) / block_size;
        let cfg = LaunchConfig {
            block_dim: (block_size, 1, 1),
            grid_dim: (grid_size, 1, 1),
            shared_mem_bytes: 0,
        };
        let n_rows_i = n_rows as i32;
        let n_cols_i = n_cols as i32;
        let mut builder = gpu.stream.launch_builder(&gpu.compute_row_norms_fn);
        builder.arg(&d_data);
        builder.arg(&mut d_norms);
        builder.arg(&n_rows_i);
        builder.arg(&n_cols_i);
        unsafe { builder.launch(cfg) }.context("Failed to launch compute_row_norms")?;
    }
    let norms: Vec<f32> = gpu.stream.clone_dtoh(&d_norms)?;

    if normalize {
        // Normalize rows in-place (classic cosine path)
        let block_size = 256u32;
        let grid_y = ((n_cols as u32) + block_size - 1) / block_size;
        let cfg = LaunchConfig {
            block_dim: (block_size, 1, 1),
            grid_dim: (n_rows as u32, grid_y, 1),
            shared_mem_bytes: 0,
        };
        let n_rows_i = n_rows as i32;
        let n_cols_i = n_cols as i32;
        let mut builder = gpu.stream.launch_builder(&gpu.row_l2_normalize_fn);
        builder.arg(&mut d_data);
        builder.arg(&d_norms);
        builder.arg(&n_rows_i);
        builder.arg(&n_cols_i);
        unsafe { builder.launch(cfg) }.context("Failed to launch row_l2_normalize")?;
    }

    // Compute similarity matrix S = A * A^T via cuBLAS sgemm
    // Data is row-major (n_rows x n_cols). cuBLAS expects column-major.
    // Row-major A is col-major A^T (n_cols x n_rows).
    // We want C = A * A^T. In col-major terms with our stored matrix B = A^T (n_cols x n_rows):
    // C = B^T * B, so we call sgemm with transa=T, transb=N on B.
    // C is n_rows x n_rows, stored column-major.
    // m = n_rows (rows of op(A)), n = n_rows (cols of op(B)), k = n_cols (cols of op(A))
    // lda = n_cols (leading dim of B before transpose), ldb = n_cols, ldc = n_rows
    let mut d_sim = gpu
        .stream
        .alloc_zeros::<f32>(n_rows * n_rows)
        .context("Failed to allocate similarity matrix")?;

    unsafe {
        gpu.cublas.gemm(
            GemmConfig {
                transa: sys::cublasOperation_t::CUBLAS_OP_T,
                transb: sys::cublasOperation_t::CUBLAS_OP_N,
                m: n_rows as i32,
                n: n_rows as i32,
                k: n_cols as i32,
                alpha: 1.0f32,
                lda: n_cols as i32,
                ldb: n_cols as i32,
                beta: 0.0f32,
                ldc: n_rows as i32,
            },
            &d_data,  // A (col-major storage = row-major A)
            &d_data,  // B (same)
            &mut d_sim,
        )
    }
    .context("cuBLAS sgemm failed")?;

    // Download similarity matrix to CPU
    log::info!("Downloading similarity matrix from GPU...");
    let sim: Vec<f32> = gpu
        .stream
        .clone_dtoh(&d_sim)
        .context("Failed to download similarity matrix")?;

    // Extract top-k neighbors per row (excluding self)
    // sim is column-major n_rows x n_rows: element (i,j) is at sim[i + j*n_rows]
    log::info!("Extracting top-{} neighbors...", k);
    let (indices, distances) = extract_topk(&sim, n_rows, k, if normalize { None } else { Some(&norms) });

    Ok(KnnGraph {
        n_points: n_rows,
        n_neighbors: k,
        indices,
        distances,
    })
}

pub fn compute_hamming_knn(
    gpu: &GpuContext,
    data: &[f32],
    n_rows: usize,
    n_cols: usize,
    k: usize,
) -> Result<KnnGraph> {
    log::info!(
        "Computing Hamming k-NN on GPU: {} x {}, k={}",
        n_rows,
        n_cols,
        k
    );

    let d_data = gpu
        .stream
        .clone_htod(data)
        .context("Failed to upload data to GPU")?;

    let mut d_dist = gpu
        .stream
        .alloc_zeros::<f32>(n_rows * n_rows)
        .context("Failed to allocate distance matrix")?;

    {
        let total = (n_rows * n_rows) as u32;
        let block_size = 256u32;
        let grid_size = (total + block_size - 1) / block_size;
        let cfg = LaunchConfig {
            block_dim: (block_size, 1, 1),
            grid_dim: (grid_size, 1, 1),
            shared_mem_bytes: 0,
        };
        let n_rows_i = n_rows as i32;
        let n_cols_i = n_cols as i32;
        let mut builder = gpu.stream.launch_builder(&gpu.compute_hamming_matrix_fn);
        builder.arg(&d_data);
        builder.arg(&mut d_dist);
        builder.arg(&n_rows_i);
        builder.arg(&n_cols_i);
        unsafe { builder.launch(cfg) }.context("Failed to launch compute_hamming_matrix")?;
    }

    log::info!("Downloading distance matrix from GPU...");
    let dist: Vec<f32> = gpu
        .stream
        .clone_dtoh(&d_dist)
        .context("Failed to download distance matrix")?;

    log::info!("Extracting top-{} neighbors...", k);
    let (indices, distances) = extract_topk_dist(&dist, n_rows, k);

    Ok(KnnGraph {
        n_points: n_rows,
        n_neighbors: k,
        indices,
        distances,
    })
}

/// Extract top-k from a distance matrix (smallest distances = nearest neighbors)
fn extract_topk_dist(dist: &[f32], n: usize, k: usize) -> (Vec<Vec<usize>>, Vec<Vec<f32>>) {
    let (indices, distances): (Vec<_>, Vec<_>) = (0..n)
        .into_par_iter()
        .map(|i| {
            // Row-major: element (i,j) at dist[i * n + j]
            let mut pairs: Vec<(f32, usize)> = (0..n)
                .filter(|&j| j != i)
                .map(|j| (dist[i * n + j], j))
                .collect();

            // Partial sort: smallest distances first
            if pairs.len() > k {
                pairs.select_nth_unstable_by(k - 1, |a, b| {
                    a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
                });
                pairs.truncate(k);
            }
            pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

            let idx: Vec<usize> = pairs.iter().map(|p| p.1).collect();
            let d: Vec<f32> = pairs.iter().map(|p| p.0.max(0.0)).collect();
            (idx, d)
        })
        .unzip();

    (indices, distances)
}

fn extract_topk(sim: &[f32], n: usize, k: usize, norms: Option<&[f32]>) -> (Vec<Vec<usize>>, Vec<Vec<f32>>) {
    let (indices, distances): (Vec<_>, Vec<_>) = (0..n)
        .into_par_iter()
        .map(|i| {
            // Collect (similarity, index) pairs, excluding self
            // Column-major: element (i,j) at sim[i + j*n]
            // If norms provided, normalize dot products to cosine: sim / (norm_i * norm_j)
            let norm_i = norms.map(|n| n[i]).unwrap_or(1.0);
            let mut pairs: Vec<(f32, usize)> = (0..n)
                .filter(|&j| j != i)
                .map(|j| {
                    let dot = sim[i + j * n];
                    let cosine = match norms {
                        Some(ns) => {
                            let denom = norm_i * ns[j];
                            if denom > 1e-10 { dot / denom } else { 0.0 }
                        }
                        None => dot,
                    };
                    (cosine, j)
                })
                .collect();

            // Partial sort to get top-k by similarity (highest first)
            if pairs.len() > k {
                pairs.select_nth_unstable_by(k - 1, |a, b| {
                    b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
                });
                pairs.truncate(k);
            }
            pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

            let idx: Vec<usize> = pairs.iter().map(|p| p.1).collect();
            // Convert cosine similarity to cosine distance
            let dist: Vec<f32> = pairs.iter().map(|p| (1.0 - p.0).max(0.0)).collect();
            (idx, dist)
        })
        .unzip();

    (indices, distances)
}
