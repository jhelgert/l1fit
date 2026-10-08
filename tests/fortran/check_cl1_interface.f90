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
