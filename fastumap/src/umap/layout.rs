use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;

use super::fuzzy_set::FuzzyGraph;

pub struct LayoutParams {
    pub n_epochs: usize,
    pub learning_rate: f32,
    pub min_dist: f32,
    pub spread: f32,
    pub negative_sample_rate: usize,
}

pub fn optimize_layout(
    embedding: &mut Vec<[f32; 2]>,
    graph: &FuzzyGraph,
    params: &LayoutParams,
    seed: u64,
) {
    let n = graph.n_points;
    let (a, b) = find_ab_params(params.spread, params.min_dist);
    log::info!(
        "Optimizing layout: {} epochs, a={:.4}, b={:.4}",
        params.n_epochs,
        a,
        b
    );

    // Compute epochs_per_sample for each edge
    let max_weight = graph
        .edges
        .iter()
        .map(|e| e.weight)
        .fold(0.0f32, f32::max);

    // epochs_per_sample: how many epochs between samples of this edge.
    // Highest-weight edges are sampled every epoch (=1.0).
    // Lower-weight edges are sampled less frequently.
    let epochs_per_sample: Vec<f32> = graph
        .edges
        .iter()
        .map(|e| {
            if e.weight > 0.0 {
                max_weight / e.weight
            } else {
                f32::INFINITY
            }
        })
        .collect();

    let epochs_per_negative_sample: Vec<f32> = epochs_per_sample
        .iter()
        .map(|&e| e / params.negative_sample_rate as f32)
        .collect();

    let mut epoch_of_next_sample: Vec<f32> = epochs_per_sample.clone();
    let mut epoch_of_next_negative_sample: Vec<f32> = epochs_per_negative_sample.clone();

    let mut rng = Pcg64::seed_from_u64(seed);

    for epoch in 0..params.n_epochs {
        let alpha = params.learning_rate * (1.0 - epoch as f32 / params.n_epochs as f32);

        for (edge_idx, edge) in graph.edges.iter().enumerate() {
            if epoch_of_next_sample[edge_idx] > epoch as f32 {
                continue;
            }

            let i = edge.row;
            let j = edge.col;

            let dx = embedding[i][0] - embedding[j][0];
            let dy = embedding[i][1] - embedding[j][1];
            let dist_sq = dx * dx + dy * dy;

            // Attractive force: gradient of UMAP kernel w.r.t. embedding coords
            let grad_coeff = if dist_sq > 0.0 {
                -2.0 * a * b * dist_sq.powf(b - 1.0) / (1.0 + a * dist_sq.powf(b))
            } else {
                0.0
            };
            let grad_x = clip(grad_coeff * dx);
            let grad_y = clip(grad_coeff * dy);

            embedding[i][0] += alpha * grad_x;
            embedding[i][1] += alpha * grad_y;
            embedding[j][0] -= alpha * grad_x;
            embedding[j][1] -= alpha * grad_y;

            // Repulsive forces (negative sampling)
            let n_neg = ((epoch as f32 - epoch_of_next_negative_sample[edge_idx])
                / epochs_per_negative_sample[edge_idx])
                .floor() as usize;
            let n_neg = n_neg.min(params.negative_sample_rate * 2);

            for _ in 0..n_neg {
                let k = rng.gen_range(0..n);
                if k == i {
                    continue;
                }

                let dx = embedding[i][0] - embedding[k][0];
                let dy = embedding[i][1] - embedding[k][1];
                let dist_sq = dx * dx + dy * dy;

                let grad_coeff = if dist_sq > 0.0 {
                    2.0 * b / ((0.001 + dist_sq) * (1.0 + a * dist_sq.powf(b)))
                } else {
                    0.0
                };
                let grad_x = clip(grad_coeff * dx);
                let grad_y = clip(grad_coeff * dy);

                embedding[i][0] += alpha * grad_x;
                embedding[i][1] += alpha * grad_y;
            }

            epoch_of_next_sample[edge_idx] += epochs_per_sample[edge_idx];
            epoch_of_next_negative_sample[edge_idx] += epochs_per_negative_sample[edge_idx];
        }

        if epoch % 50 == 0 || epoch == params.n_epochs - 1 {
            log::info!("Epoch {}/{}", epoch + 1, params.n_epochs);
        }
    }
}

fn clip(val: f32) -> f32 {
    val.clamp(-4.0, 4.0)
}

fn find_ab_params(spread: f32, min_dist: f32) -> (f32, f32) {
    // Curve fit: (1 + a * d^(2b))^(-1) should approximate
    // 1.0 for d <= min_dist, exp(-(d - min_dist) / spread) for d > min_dist
    //
    // Use f64 for fitting precision, then convert back to f32.
    let spread = spread as f64;
    let min_dist = min_dist as f64;

    let xs: Vec<f64> = (1..300).map(|i| i as f64 * 0.01).collect();
    let targets: Vec<f64> = xs
        .iter()
        .map(|&d| {
            if d <= min_dist {
                1.0
            } else {
                (-(d - min_dist) / spread).exp()
            }
        })
        .collect();

    // Fine grid search
    let mut best_a = 1.0f64;
    let mut best_b = 1.0f64;
    let mut best_err = f64::INFINITY;

    for ai in 1..500 {
        let a = ai as f64 * 0.05;
        for bi in 1..200 {
            let b = bi as f64 * 0.01;
            let err: f64 = xs
                .iter()
                .zip(&targets)
                .map(|(&d, &t)| {
                    let y = 1.0 / (1.0 + a * d.powf(2.0 * b));
                    (y - t) * (y - t)
                })
                .sum();
            if err < best_err {
                best_err = err;
                best_a = a;
                best_b = b;
            }
        }
    }

    log::info!(
        "find_ab_params(spread={}, min_dist={}): a={:.4}, b={:.4}, err={:.6}",
        spread,
        min_dist,
        best_a,
        best_b,
        best_err
    );

    (best_a as f32, best_b as f32)
}
