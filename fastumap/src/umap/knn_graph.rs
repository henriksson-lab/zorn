pub struct KnnGraph {
    pub n_points: usize,
    pub n_neighbors: usize,
    pub indices: Vec<Vec<usize>>,
    pub distances: Vec<Vec<f32>>,
}
