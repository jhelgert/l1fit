# GPU notes (OpenACC / CUDA Fortran)

Only apply when the user asks for GPU support or the legacy code already has directives.
Modernize the serial code **first**, verify it, then add directives as a separate step.

## OpenACC
Directives are comments on non-GPU compilers, so the code stays portable.

```fortran
!$acc parallel loop collapse(2) default(present)
do j = 1, ny
   do i = 1, nx
      flux(i,j) = hu(i,j)**2 / h(i,j) + 0.5_wp * GRAVITY * h(i,j)**2
   end do
end do
```

- `parallel loop` parallelizes; `collapse(n)` merges nested loops; `default(present)` assumes data is on device; `data copyin/copy/copyout` moves data.
- **Reductions**: `!$acc parallel loop reduction(+:total)`.
- **Coalescing**: Fortran is column-major - the innermost loop must run over the *first* index (`do j ... do i ... a(i,j)`). Interchange legacy loops that stride the wrong way only if it does not change results (it can change summation order).

## Data management
Copy to the device once at solver start, keep there, copy back once at the end:

```fortran
!$acc enter data copyin(state%h, state%hu, state%hv, state%z)
do while (state%time < config%t_end)
   call compute_timestep(state, config)   ! all device work
end do
!$acc exit data copyout(state%h, state%hu, state%hv)
!$acc exit data delete(state%z)
```

Unit tests wrap kernel calls in explicit `!$acc data copyin(...) copyout(...)` regions and verify on the host.
GPU memory needs explicit cleanup (`exit data`) - types with device data are the exception to "no destroy routines".

## CUDA Fortran
```fortran
attributes(global) subroutine compute_flux_kernel(h, hu, flux, nx, ny)
   real(wp), device, intent(in)  :: h(nx, ny), hu(nx, ny)
   real(wp), device, intent(out) :: flux(nx, ny)
   integer, value, intent(in) :: nx, ny
   integer :: i, j
   i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
   j = (blockIdx%y - 1) * blockDim%y + threadIdx%y
   if (i <= nx .and. j <= ny) flux(i,j) = hu(i,j)**2 / h(i,j) + 0.5_wp * GRAVITY * h(i,j)**2
end subroutine
```
- `attributes(global)` = kernel callable from host; `attributes(device)` = device-only; `device` = GPU memory; scalars passed to kernels need `value`.
- Reductions: `atomicAdd` for simple cases, parallel reduction kernels for performance-critical code.

## Multi-GPU with MPI
One rank per GPU, bound by node-local rank:

```fortran
call MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_comm, ierr)
call MPI_Comm_rank(node_comm, node_rank, ierr)
!$acc set device_num(node_rank)        ! or ierr = cudaSetDevice(node_rank)
```
