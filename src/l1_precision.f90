module l1_precision
   !! Working precision of the L1 fitting code.
   !!
   !! `wp` is the working precision of the arrays and scalars of the solver; the original
   !! ACM TOMS algorithm 552 used single precision `REAL` here. `dp` is the precision of the
   !! accumulator that sums the residuals (the original `DOUBLE PRECISION`).
   !! The kinds are taken from `iso_c_binding` because the solver is called from C / Python
   !! through a `bind(C)` entry point.
   !! Changing `wp` here changes the precision of the whole library.
   use, intrinsic :: iso_c_binding, only: c_double
   implicit none (type, external)
   private

   public :: wp, dp

   integer, parameter :: wp = c_double
      !! Working precision (IEEE binary64; the original code used `c_float`, i.e. binary32)
   integer, parameter :: dp = c_double
      !! Accumulator precision (kind of the original `DOUBLE PRECISION`, IEEE binary64)

end module l1_precision
