use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use rayon::prelude::*;
use std::collections::HashMap;

use super::knn_graph::KnnGraph;

pub struct SparseEdge {
    pub row: usize,
    pub col: usize,
    pub weight: f32,
}

pub struct FuzzyGraph {
    pub n_points: usize,
    pub edges: Vec<SparseEdge>,
}

pub fn fuzzy_simplicial_set(knn: &KnnGraph, seed: u64) -> FuzzyGraph {
    let n = knn.n_points;
    let k = knn.n_neighbors;
    let target = (k as f32).ln() / (2.0f32).ln(); // log2(k)

    log::info!(
        "Building fuzzy simplicial set: {} points, target={}",
        n,
        target
    );

    // Compute rho (distance to nearest neighbor) and sigma per point
    let (rhos, sigmas): (Vec<f32>, Vec<f32>) = (0..n)
        .into_par_iter()
        .map(|i| {
            let dists = &knn.distances[i];
            let rho = dists.iter().cloned().fold(f32::INFINITY, f32::min);
            let rho = if rho.is_finite() { rho } else { 0.0 };
            let sigma = find_sigma(dists, rho, target);
            (rho, sigma)
        })
        .unzip();

    // Compute directed membership strengths
    let mut edge_map: HashMap<(usize, usize), f32> = HashMap::new();

    for i in 0..n {
        for jj in 0..k {
            let j = knn.indices[i][jj];
            let d = knn.distances[i][jj];
            let w = if d <= rhos[i] {
                1.0
            } else if sigmas[i] > 1e-10 {
                (-(d - rhos[i]) / sigmas[i]).exp()
            } else {
                0.0
            };
            if w > 1e-10 {
                edge_map.insert((i, j), w);
            }
        }
    }

    // Symmetrize: w_sym = w_ij + w_ji - w_ij * w_ji
    let mut sym_map: HashMap<(usize, usize), f32> = HashMap::new();
    for (&(i, j), &w_ij) in &edge_map {
        let w_ji = edge_map.get(&(j, i)).copied().unwrap_or(0.0);
        let w_sym = w_ij + w_ji - w_ij * w_ji;
        if w_sym > 1e-10 {
            let key = if i <= j { (i, j) } else { (j, i) };
            sym_map.insert(key, w_sym);
        }
    }

    // Build both directions, then sort by (row, col) to match reference COO ordering
    let mut edges: Vec<SparseEdge> = Vec::with_capacity(sym_map.len() * 2);
    for ((i, j), w) in &sym_map {
        edges.push(SparseEdge {
            row: *i,
            col: *j,
            weight: *w,
        });
        edges.push(SparseEdge {
            row: *j,
            col: *i,
            weight: *w,
        });
    }
    // Deterministic shuffle to distribute updates evenly across points
    let mut rng = Pcg64::seed_from_u64(seed + 100);
    edges.shuffle(&mut rng);

    log::info!("Fuzzy graph: {} edges", edges.len());

    FuzzyGraph {
        n_points: n,
        edges,
    }
}

fn find_sigma(distances: &[f32], rho: f32, target: f32) -> f32 {
    let mut lo = 1e-10f32;
    let mut hi = 1000.0f32;
    let mut mid = 1.0f32;

    for _ in 0..64 {
        mid = (lo + hi) / 2.0;
        let val: f32 = distances
            .iter()
            .map(|&d| {
                let shifted = d - rho;
                if shifted <= 0.0 {
                    1.0
                } else {
                    (-shifted / mid).exp()
                }
            })
            .sum();

        if (val - target).abs() < 1e-5 {
            break;
        }
        if val > target {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    mid
}
