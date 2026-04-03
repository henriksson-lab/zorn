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
extern "C" __global__ void compute_hamming_matrix(
    const float* data,
    float* dist,
    int n_rows,
    int n_cols
) {
    // Each thread computes one (i, j) pair
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n_rows * n_rows;
    if (idx >= total) return;

    int i = idx / n_rows;
    int j = idx % n_rows;

    if (i == j) {
        dist[idx] = 0.0f;
        return;
    }

    int mismatches = 0;
    int active = 0;
    for (int c = 0; c < n_cols; c++) {
        float a = data[i * n_cols + c];
        float b = data[j * n_cols + c];
        int a_nz = (a != 0.0f);
        int b_nz = (b != 0.0f);
        if (a_nz || b_nz) {
            active++;
            // sign mismatch: different signs, or one is zero and other isn't
            float sa = (a > 0.0f) ? 1.0f : ((a < 0.0f) ? -1.0f : 0.0f);
            float sb = (b > 0.0f) ? 1.0f : ((b < 0.0f) ? -1.0f : 0.0f);
            if (sa != sb) mismatches++;
        }
    }

    dist[idx] = (active > 0) ? ((float)mismatches / (float)active) : 1.0f;
}
"#;
