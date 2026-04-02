pub const KERNELS_SRC: &str = r#"
extern "C" __global__ void compute_row_norms(
    const float* data,
    float* norms,
    int n_rows,
    int n_cols
) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n_rows) return;

    float sum = 0.0f;
    for (int col = 0; col < n_cols; col++) {
        float val = data[row * n_cols + col];
        sum += val * val;
    }
    norms[row] = sqrtf(sum);
}

extern "C" __global__ void row_l2_normalize(
    float* data,
    const float* norms,
    int n_rows,
    int n_cols
) {
    int row = blockIdx.x;
    int col = threadIdx.x + blockIdx.y * blockDim.x;
    if (row >= n_rows || col >= n_cols) return;

    float norm = norms[row];
    if (norm > 1e-10f) {
        data[row * n_cols + col] /= norm;
    }
}
"#;
