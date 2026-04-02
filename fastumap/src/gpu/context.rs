use anyhow::{Context, Result};
use cudarc::cublas::safe::CudaBlas;
use cudarc::driver::{CudaContext, CudaFunction, CudaModule, CudaStream};
use std::sync::Arc;

use super::kernels::KERNELS_SRC;

pub struct GpuContext {
    pub ctx: Arc<CudaContext>,
    pub stream: Arc<CudaStream>,
    pub cublas: CudaBlas,
    pub module: Arc<CudaModule>,
    pub compute_row_norms_fn: CudaFunction,
    pub row_l2_normalize_fn: CudaFunction,
}

impl GpuContext {
    pub fn new(device_id: usize) -> Result<Self> {
        let ctx = CudaContext::new(device_id).context("Failed to create CUDA context")?;
        let stream = ctx.default_stream();

        let cublas =
            CudaBlas::new(stream.clone()).context("Failed to create cuBLAS handle")?;

        let ptx = cudarc::nvrtc::compile_ptx(KERNELS_SRC)
            .context("Failed to compile CUDA kernels")?;

        let module = ctx.load_module(ptx).context("Failed to load CUDA module")?;

        let compute_row_norms_fn = module
            .load_function("compute_row_norms")
            .context("Failed to load compute_row_norms kernel")?;
        let row_l2_normalize_fn = module
            .load_function("row_l2_normalize")
            .context("Failed to load row_l2_normalize kernel")?;

        Ok(Self {
            ctx,
            stream,
            cublas,
            module,
            compute_row_norms_fn,
            row_l2_normalize_fn,
        })
    }
}
