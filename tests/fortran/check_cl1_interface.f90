program check_cl1_interface
   !! Checks of the Fortran interface `cl1` that cannot be reached through the C ABI, because the C
   !! wrapper derives the array shapes from the dimension arguments: argument validation and the
   !! acceptance of arrays that are larger than required.
   !!
   !! Exit status 0 if every check passes; a message is printed and the status is 1 otherwise.
   use l1_calgo552, only: cl1, L1_OPTIMAL, L1_INVALID_INPUT
   use l1_precision, only: wp
   implicit none (type, external)

   integer, parameter :: NUM_ROWS = 5, NUM_COLS = 2
      !! Problem used throughout: fit a line `b = c0 + c1 t` to five points (k = 5, n = 2)
   integer :: num_failures

   num_failures = 0

   call check_valid_problem_is_solved()
   call check_oversized_arrays_give_identical_results()
   call check_outputs_are_written_where_documented()
   call check_every_residual_entry_is_zeroed()
   call check_invalid_arguments_are_rejected()

   if (num_failures /= 0) then
      print "(a,i0,a)", "FAILED: ", num_failures, " check(s)"
      stop 1
   end if
   print "(a)", "all interface checks passed"

contains

   subroutine report(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print "(2a)", "check failed: ", message
         num_failures = num_failures + 1
      end if
   end subroutine report

   subroutine fill_line_fit(q)
      !! Five points (t, b) = (0, 1), (1, 3), (2, 5), (3, 7), (4, 30): the last one is an outlier.
      real(wp), intent(out) :: q(:, :)
      integer :: i

      q = 0.0_wp
      do i = 1, NUM_ROWS
         q(i, 1) = 1.0_wp
         q(i, 2) = real(i - 1, wp)
      end do
      q(1:NUM_ROWS, NUM_COLS + 1) = [1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 30.0_wp]
   end subroutine fill_line_fit

   subroutine solve_exact(x, error, kode)
      real(wp), intent(out) :: x(NUM_COLS), error
      integer, intent(out) :: kode
      real(wp) :: q(NUM_ROWS + 2, NUM_COLS + 2), res(NUM_ROWS)
      integer :: iter

      call fill_line_fit(q)
      kode = 0
      iter = 100
      x = 0.0_wp
      res = 0.0_wp
      call cl1(NUM_ROWS, 0, 0, NUM_COLS, q, kode, 1.0e-10_wp, iter, x, res, error)
   end subroutine solve_exact

   subroutine check_valid_problem_is_solved()
      real(wp) :: x(NUM_COLS), error
      integer :: kode

      call solve_exact(x, error, kode)
      call report(kode == L1_OPTIMAL, "valid problem: kode == L1_OPTIMAL")
      ! The L1 fit ignores the outlier: the line through the first four points, error = 30 - 9 = 21
      call report(all(abs(x - [1.0_wp, 2.0_wp]) < 1.0e-9_wp), "valid problem: x == (1, 2)")
      call report(abs(error - 21.0_wp) < 1.0e-9_wp, "valid problem: error == 21")
   end subroutine check_valid_problem_is_solved

   subroutine check_oversized_arrays_give_identical_results()
      !! Arrays may be larger than required (legacy contract); results must be bit-identical.
      real(wp) :: x_exact(NUM_COLS), error_exact
      real(wp) :: q(NUM_ROWS + 9, NUM_COLS + 7), x(NUM_COLS + 5), res(NUM_ROWS + 3), error
      integer :: kode_exact, kode, iter

      call solve_exact(x_exact, error_exact, kode_exact)

      q = 0.0_wp
      call fill_line_fit(q(1:NUM_ROWS + 2, 1:NUM_COLS + 2))
      x = 0.0_wp
      res = 0.0_wp
      kode = 0
      iter = 100
      call cl1(NUM_ROWS, 0, 0, NUM_COLS, q, kode, 1.0e-10_wp, iter, x, res, error)

      call report(kode == kode_exact, "oversized arrays: same kode")
      call report(all(x(1:NUM_COLS) == x_exact), "oversized arrays: bit-identical x")
      call report(error == error_exact, "oversized arrays: bit-identical error")
   end subroutine check_oversized_arrays_give_identical_results

   subroutine check_outputs_are_written_where_documented()
      !! `x(1:n)` and `res(1:k+l+m)` are completely overwritten, including the entries of variables
      !! and rows that are not in the final basis (these must read exactly 0, not a stale value),
      !! and nothing beyond them is touched.
      !!
      !! Which rows end up outside the basis depends on the data, so the same problem is solved with
      !! the rows cyclically shifted: every row position gets a turn at not being in the basis, and a
      !! missing zeroing of any single `res` entry shows up as a stale sentinel value.
      real(wp), parameter :: SENTINEL = 12345.0_wp
      integer, parameter :: NUM_VARS = 3
      real(wp), parameter :: T(NUM_ROWS) = [4.0_wp, 0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp]
      real(wp), parameter :: B(NUM_ROWS) = [30.0_wp, 1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp]
      real(wp) :: q(NUM_ROWS + 2, NUM_VARS + 2), x(NUM_VARS + 4), res(NUM_ROWS + 4), error
      real(wp) :: expected_res(NUM_ROWS)
      integer :: kode, iter, shift
      character(len=32) :: label

      do shift = 0, NUM_ROWS - 1
         write (label, "(a,i0)") "outputs, rows shifted by ", shift
         q = 0.0_wp
         q(1:NUM_ROWS, 1) = 1.0_wp
         q(1:NUM_ROWS, 2) = cshift(T, shift)
         ! The third column stays 0: the third variable never enters the basis.
         q(1:NUM_ROWS, NUM_VARS + 1) = cshift(B, shift)
         ! Row 1 of the unshifted data is the outlier (residual 21); all other rows fit exactly.
         expected_res = 0.0_wp
         expected_res(1 + modulo(-shift, NUM_ROWS)) = 21.0_wp

         x = SENTINEL
         res = SENTINEL
         kode = 0
         iter = 100
         call cl1(NUM_ROWS, 0, 0, NUM_VARS, q, kode, 1.0e-10_wp, iter, x, res, error)

         call report(kode == L1_OPTIMAL, trim(label)//": kode == L1_OPTIMAL")
         call report(all(abs(x(1:2) - [1.0_wp, 2.0_wp]) < 1.0e-9_wp), trim(label)//": x(1:2) == (1, 2)")
         call report(x(3) == 0.0_wp, trim(label)//": the variable outside the basis reads exactly 0")
         call report(all(x(NUM_VARS + 1:) == SENTINEL), trim(label)//": x beyond n is untouched")
         call report(all(abs(res(1:NUM_ROWS) - expected_res) < 1.0e-9_wp), trim(label)//": res(1:k) == b - A x")
         call report(all(res(NUM_ROWS + 1:) == SENTINEL), trim(label)//": res beyond k+l+m is untouched")
      end do
   end subroutine check_outputs_are_written_where_documented

   subroutine check_every_residual_entry_is_zeroed()
      !! The L1 fit of one constant to three numbers is their median. The row holding the median is
      !! interpolated exactly, so its residual variable is not in the final basis and is never
      !! assigned by the solver: it reads 0 only because the outputs are zeroed first. Putting the
      !! median in each of the three row positions makes every `res` entry unassigned once.
      real(wp), parameter :: SENTINEL = 12345.0_wp
      integer, parameter :: NUM_POINTS = 3
      real(wp), parameter :: VALUES(NUM_POINTS, NUM_POINTS) = reshape( &
         [2.0_wp, 1.0_wp, 3.0_wp, &   ! median in row 1
          1.0_wp, 2.0_wp, 3.0_wp, &   ! median in row 2
          1.0_wp, 3.0_wp, 2.0_wp], [NUM_POINTS, NUM_POINTS])   ! median in row 3
      real(wp) :: q(NUM_POINTS + 2, 3), x(1), res(NUM_POINTS), error
      integer :: kode, iter, median_row

      do median_row = 1, NUM_POINTS
         q = 0.0_wp
         q(1:NUM_POINTS, 1) = 1.0_wp
         q(1:NUM_POINTS, 2) = VALUES(:, median_row)
         x = SENTINEL
         res = SENTINEL
         kode = 0
         iter = 100
         call cl1(NUM_POINTS, 0, 0, 1, q, kode, 1.0e-10_wp, iter, x, res, error)

         call report(kode == L1_OPTIMAL, "median: kode == L1_OPTIMAL")
         call report(abs(x(1) - 2.0_wp) < 1.0e-9_wp, "median: x == 2")
         call report(res(median_row) == 0.0_wp, "median: the interpolated row reads exactly 0")
         call report(all(abs(res - (VALUES(:, median_row) - 2.0_wp)) < 1.0e-9_wp), "median: res == b - x")
      end do
   end subroutine check_every_residual_entry_is_zeroed

   subroutine check_rejected(k, l, m, n, q, x, res, label)
      !! Call `cl1` with arguments that must be rejected and check the reported state.
      integer, intent(in) :: k, l, m, n
      real(wp), contiguous, intent(inout) :: q(:, :), x(:), res(:)
      character(len=*), intent(in) :: label
      real(wp) :: error
      integer :: kode, iter

      kode = 0
      iter = 7
      error = -1.0_wp
      call cl1(k, l, m, n, q, kode, 1.0e-10_wp, iter, x, res, error)
      call report(kode == L1_INVALID_INPUT, "invalid input is rejected: "//label)
      call report(iter == 0 .and. error == 0.0_wp, "invalid input: iter and error reset: "//label)
   end subroutine check_rejected

   subroutine check_invalid_arguments_are_rejected()
      !! Every case violates exactly one condition: all other arrays are generously sized, so that
      !! a weakened (or missing) check for that condition cannot be masked by another one.
      integer, parameter :: BIG = 20
      real(wp) :: q(BIG, BIG), x(BIG), res(BIG), q_before(BIG, BIG)
      real(wp) :: q_short_k(NUM_ROWS + 1, BIG), q_short_klm(NUM_ROWS + 3, BIG)
         !! contiguous arrays that are exactly one row (k) or two rows (k, l, m counted) too short

      q = 0.0_wp
      q_short_k = 0.0_wp
      q_short_klm = 0.0_wp
      call fill_line_fit(q(1:NUM_ROWS + 2, 1:NUM_COLS + 2))
      q_before = q
      x = 0.0_wp
      res = 0.0_wp

      ! Dimensions
      call check_rejected(0, 0, 0, NUM_COLS, q, x, res, "k = 0")
      call check_rejected(NUM_ROWS, -1, 0, NUM_COLS, q, x, res, "l < 0")
      call check_rejected(NUM_ROWS, 0, -1, NUM_COLS, q, x, res, "m < 0")
      call check_rejected(NUM_ROWS, 0, 0, 0, q, x, res, "n = 0")
      ! Sizes, with k = 5, n = 2 (q needs 7 x 4, x needs 2, res needs 5) ...
      call check_rejected(NUM_ROWS, 0, 0, NUM_COLS, q_short_k, x, res, "q one row short")
      call check_rejected(NUM_ROWS, 0, 0, NUM_COLS, q(:, 1:NUM_COLS + 1), x, res, "q one column short")
      call check_rejected(NUM_ROWS, 0, 0, NUM_COLS, q, x(1:NUM_COLS - 1), res, "x one element short")
      call check_rejected(NUM_ROWS, 0, 0, NUM_COLS, q, x, res(1:NUM_ROWS - 1), "res one element short")
      ! ... and with all of k, l and m counted (q needs 9 rows, res needs 7 elements)
      call check_rejected(NUM_ROWS, 1, 1, NUM_COLS, q_short_klm, x, res, "q short by m and l")
      call check_rejected(NUM_ROWS, 1, 1, NUM_COLS, q, x, res(1:NUM_ROWS + 1), "res short by m")
      call check_rejected(NUM_ROWS, 1, 0, NUM_COLS, q, x, res(1:NUM_ROWS), "res short by l")

      call report(all(q == q_before), "invalid input: q is left untouched")
   end subroutine check_invalid_arguments_are_rejected

end program check_cl1_interface
