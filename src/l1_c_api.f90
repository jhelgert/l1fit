module l1_c_api
   !! C-callable entry point of the L1 solver (used by the Python / C extension and by `ctypes`).
   !!
   !! C has no array shape information, so this thin wrapper receives the problem dimensions and
   !! declares the arrays as explicit-shape (contiguous, column-major) with exactly the sizes
   !! implied by them, then calls the Fortran interface `cl1`. Passing explicit-shape arrays on to
   !! the assumed-shape dummies of `cl1` does not copy them.
   use, intrinsic :: iso_c_binding, only: c_int
   use l1_calgo552, only: cl1
   use l1_precision, only: wp
   implicit none (type, external)
   private

   public :: l1_cl1

contains

   ! allow(too-many-arguments)
   subroutine l1_cl1(k, l, m, n, q, kode, toler, iter, x, res, error) bind(C, name="l1_cl1")
      !! C interface of `cl1`; see `cl1` for the meaning of every argument.
      !!
      !! Contract for the caller: `q` is a contiguous column-major array with `k+l+m+2` rows and
      !! `n+2` columns, `x` has `n` elements and `res` has `k+l+m` elements.
      integer(c_int), intent(in) :: k
         !! Number of rows of `A`
      integer(c_int), intent(in) :: l
         !! Number of rows of `C`
      integer(c_int), intent(in) :: m
         !! Number of rows of `E`
      integer(c_int), intent(in) :: n
         !! Number of columns of `A`, `C`, `E`
      real(wp), intent(inout) :: q(k + l + m + 2, n + 2)
         !! Problem data, destroyed by the solver
      integer(c_int), intent(inout) :: kode
         !! Entry flag (0 or 1) and exit status
      real(wp), intent(in) :: toler
         !! Tolerance
      integer(c_int), intent(inout) :: iter
         !! Maximum number of iterations on entry, number performed on exit
      real(wp), intent(inout) :: x(n)
         !! Solution (and sign restrictions on entry if `kode = 1`)
      real(wp), intent(inout) :: res(k + l + m)
         !! Residuals (and sign restrictions on entry if `kode = 1`)
      real(wp), intent(out) :: error
         !! Minimum sum of absolute values of the residuals `b - A x`

      call cl1(k, l, m, n, q, kode, toler, iter, x, res, error)
   end subroutine l1_cl1

end module l1_c_api
