use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;

use super::fuzzy_set::FuzzyGraph;

pub fn spectral_layout(graph: &FuzzyGraph, seed: u64) -> Vec<[f32; 2]> {
    let n = graph.n_points;
    log::info!("Computing spectral initialization for {} points...", n);

    // Build adjacency: for each node, collect (neighbor, weight)
    let mut adj: Vec<Vec<(usize, f32)>> = vec![Vec::new(); n];
    for e in &graph.edges {
        adj[e.row].push((e.col, e.weight));
    }

    // Compute degree
    let degree: Vec<f32> = adj
        .iter()
        .map(|neighbors| neighbors.iter().map(|(_, w)| w).sum::<f32>())
        .collect();

    let d_inv_sqrt: Vec<f32> = degree
        .iter()
        .map(|&d| if d > 1e-10 { 1.0 / d.sqrt() } else { 0.0 })
        .collect();

    // Power iteration for eigenvectors of D^{-1/2} * W * D^{-1/2}
    // We need eigenvectors 2 and 3 (skipping the trivial constant eigenvector #1)
    let mut rng = Pcg64::seed_from_u64(seed);

    // Find 3 eigenvectors, then discard the first (trivial)
    let n_vecs = 3;
    let mut vecs: Vec<Vec<f32>> = (0..n_vecs)
        .map(|_| (0..n).map(|_| rng.gen::<f32>() - 0.5).collect())
        .collect();

    for _iter in 0..500 {
        for vi in 0..n_vecs {
            // y = D^{-1/2} W D^{-1/2} x
            let x = &vecs[vi];
            let mut y = vec![0.0f32; n];

            // z = D^{-1/2} * x
            let z: Vec<f32> = x.iter().zip(&d_inv_sqrt).map(|(xi, di)| xi * di).collect();

            // w = W * z
            for i in 0..n {
                let mut sum = 0.0f32;
                for &(j, wt) in &adj[i] {
                    sum += wt * z[j];
                }
                y[i] = sum;
            }

            // y = D^{-1/2} * w
            for i in 0..n {
                y[i] *= d_inv_sqrt[i];
            }

            // Gram-Schmidt orthogonalization against previous eigenvectors
            for vj in 0..vi {
                let dot: f32 = y.iter().zip(&vecs[vj]).map(|(a, b)| a * b).sum();
                for i in 0..n {
                    y[i] -= dot * vecs[vj][i];
                }
            }

            // Normalize
            let norm: f32 = y.iter().map(|v| v * v).sum::<f32>().sqrt();
            if norm > 1e-10 {
                for v in y.iter_mut() {
                    *v /= norm;
                }
            }

            vecs[vi] = y;
        }
    }

    // Use eigenvectors 1 and 2 (0-indexed), skipping eigenvector 0 (trivial)
    let ev1 = &vecs[1];
    let ev2 = &vecs[2];

    // Scale to have small spread for UMAP initialization
    let std1 = (ev1.iter().map(|v| v * v).sum::<f32>() / n as f32).sqrt();
    let std2 = (ev2.iter().map(|v| v * v).sum::<f32>() / n as f32).sqrt();
    let target_std = 1e-4;

    let scale1 = if std1 > 1e-10 { target_std / std1 } else { target_std };
    let scale2 = if std2 > 1e-10 { target_std / std2 } else { target_std };

    let mut result: Vec<[f32; 2]> = Vec::with_capacity(n);
    for i in 0..n {
        result.push([ev1[i] * scale1, ev2[i] * scale2]);
    }

    // Check for degenerate result
    let has_valid = result
        .iter()
        .any(|p| p[0].is_finite() && p[0].abs() > 1e-12);
    if !has_valid {
        log::warn!("Spectral initialization failed, falling back to random");
        return random_layout(n, seed);
    }

    result
}

pub fn random_layout(n: usize, seed: u64) -> Vec<[f32; 2]> {
    let mut rng = Pcg64::seed_from_u64(seed + 1);
    (0..n)
        .map(|_| {
            [
                (rng.gen::<f32>() - 0.5) * 20.0,
                (rng.gen::<f32>() - 0.5) * 20.0,
            ]
        })
        .collect()
}
