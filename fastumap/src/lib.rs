pub mod csv_io;
pub mod gpu;
pub mod umap;

use anyhow::Result;
use ndarray::Array2;

pub struct UmapParams {
    pub n_neighbors: usize,
    pub min_dist: f32,
    pub n_epochs: usize,
    pub learning_rate: f32,
    pub spread: f32,
    pub negative_sample_rate: usize,
    pub seed: u64,
    pub gpu_device: usize,
    pub random_init: bool,
    pub metric: String,
    pub transform: String,
}

impl Default for UmapParams {
    fn default() -> Self {
        Self {
            n_neighbors: 15,
            min_dist: 0.1,
            n_epochs: 500,
            learning_rate: 1.0,
            spread: 1.0,
            negative_sample_rate: 5,
            seed: 42,
            gpu_device: 0,
            random_init: false,
            metric: "cosine".to_string(),
            transform: "none".to_string(),
        }
    }
}

pub struct UmapResult {
    pub embedding: Vec<[f32; 2]>,
}

pub fn run_umap(data: &Array2<f32>, params: &UmapParams) -> Result<UmapResult> {
    let (n_rows, n_cols) = data.dim();
    log::info!("Running UMAP on {} x {} matrix", n_rows, n_cols);

    // Step 1: GPU-accelerated k-NN
    let gpu_ctx = gpu::context::GpuContext::new(params.gpu_device)?;
    let flat_data = data.as_slice().expect("Data must be contiguous");
    let normalize = params.transform != "sign2";
    let knn = match params.metric.as_str() {
        "hamming" => gpu::cosine_knn::compute_hamming_knn(
            &gpu_ctx,
            flat_data,
            n_rows,
            n_cols,
            params.n_neighbors,
        )?,
        _ => gpu::cosine_knn::compute_cosine_knn(
            &gpu_ctx,
            flat_data,
            n_rows,
            n_cols,
            params.n_neighbors,
            normalize,
        )?,
    };

    // Step 2: Fuzzy simplicial set
    let graph = umap::fuzzy_set::fuzzy_simplicial_set(&knn, params.seed);

    // Step 3: Initialization
    let mut embedding = if params.random_init {
        umap::spectral::random_layout(graph.n_points, params.seed)
    } else {
        umap::spectral::spectral_layout(&graph, params.seed)
    };

    // Step 4: SGD optimization
    let layout_params = umap::layout::LayoutParams {
        n_epochs: params.n_epochs,
        learning_rate: params.learning_rate,
        min_dist: params.min_dist,
        spread: params.spread,
        negative_sample_rate: params.negative_sample_rate,
    };
    umap::layout::optimize_layout(&mut embedding, &graph, &layout_params, params.seed);

    Ok(UmapResult { embedding })
}
