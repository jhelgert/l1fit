!--**--CH751--552--C:ID--1:8:1999
module l1_calgo552
   !! L1 solution of linear systems with linear constraints (ACM TOMS algorithm 552).
   !!
   !! Uses a modification of the simplex method of linear programming (I. Barrodale and
   !! F. D. K. Roberts, ACM Transactions on Mathematical Software, algorithm 552) to calculate an
   !! L1 solution of a `k` by `n` system of linear equations
   !!
   !!     A x = b
   !!
   !! subject to `l` linear equality constraints `C x = d` and `m` linear inequality constraints
   !! `E x <= f`, i.e. it solves
   !!
   !!     minimise ||A x - b||_1   subject to   C x = d,   E x <= f.
   !!
   !! All data are passed in one matrix `q` that is destroyed by the solver:
   !!
   !!     q(1:k+l+m, 1:n+1) = [ A  b ]
   !!                         [ C  d ]
   !!                         [ E  f ]
   !!
   !! The last row and the last column of `q(1:k+l+m+2, 1:n+2)` are workspace.
   use l1_precision, only: wp, dp
   implicit none (type, external)
   private

   public :: cl1
   public :: L1_OPTIMAL, L1_INFEASIBLE, L1_ROUNDING_ERRORS, L1_MAX_ITERATIONS
   public :: L1_INVALID_INPUT, L1_ALLOC_FAILED

   integer, parameter :: L1_OPTIMAL = 0
      !! Exit status: optimal solution found
   integer, parameter :: L1_INFEASIBLE = 1
      !! Exit status: no feasible solution to the constraints
   integer, parameter :: L1_ROUNDING_ERRORS = 2
      !! Exit status: calculations terminated prematurely due to rounding errors
   integer, parameter :: L1_MAX_ITERATIONS = 3
      !! Exit status: maximum number of iterations reached
   integer, parameter :: L1_INVALID_INPUT = 4
      !! Exit status: invalid dimensions or arrays that are too small (nothing was computed)
   integer, parameter :: L1_ALLOC_FAILED = 5
      !! Exit status: the internal workspace could not be allocated (nothing was computed)

contains

   ! The 11-argument interface (down from 18 in the original) is kept for now, and the routine is
   ! still one monolithic simplex loop; splitting it into procedures is planned (structure layer).
   ! allow(too-many-arguments, too-complex)
   subroutine cl1(k, l, m, n, q, kode, toler, iter, x, res, error)
      !! Solve the L1 problem described in the module documentation.
      !!
      !! Arrays may be larger than required; only the leading parts described below are used.
      !! `q`, `x` and `res` are declared `contiguous`: a non-contiguous actual argument (for
      !! example a strided section) is copied in and out by the compiler.
      !! If an argument is invalid `kode` is set to `L1_INVALID_INPUT` and no array is touched.
      integer, intent(in) :: k
         !! Number of rows of the matrix `A` (`k >= 1`)
      integer, intent(in) :: l
         !! Number of rows of the matrix `C` (`l >= 0`)
      integer, intent(in) :: m
         !! Number of rows of the matrix `E` (`m >= 0`)
      integer, intent(in) :: n
         !! Number of columns of the matrices `A`, `C`, `E` (`n >= 1`)
      real(wp), contiguous, intent(inout) :: q(:, :)
         !! Problem data on entry (layout in the module documentation, at least `k+l+m+2` rows and
         !! `n+2` columns); destroyed by the solver
      integer, intent(inout) :: kode
         !! On entry normally `0`. If set to `1`, sign restrictions are imposed implicitly instead of
         !! through explicit rows of `E x <= f`: see `x` and `res`.
         !! On exit one of `L1_OPTIMAL`, `L1_INFEASIBLE`, `L1_ROUNDING_ERRORS`, `L1_MAX_ITERATIONS`,
         !! `L1_INVALID_INPUT`, `L1_ALLOC_FAILED`
      real(wp), intent(in) :: toler
         !! Small positive tolerance: the solver cannot distinguish zero from any quantity whose
         !! magnitude does not exceed `toler`, in particular it does not pivot on such a number.
         !! Empirical evidence suggests `toler = 10**(-d*2/3)` where `d` is the number of decimal
         !! digits of accuracy available
      integer, intent(inout) :: iter
         !! On entry the maximum number of iterations allowed (suggested: `10*(k+l+m)`);
         !! on exit the number of simplex iterations performed
      real(wp), contiguous, intent(inout) :: x(:)
         !! On exit (`size(x) >= n`) a solution of the L1 problem. If `kode = 1` on entry, `x(j)`
         !! in {-1, 0, 1} on entry restricts the j-th variable to be <= 0, unrestricted or >= 0
      real(wp), contiguous, intent(inout) :: res(:)
         !! On exit (`size(res) >= k+l+m`) the residuals `b - A x` in the first `k` components,
         !! `d - C x` in the next `l` (these are 0) and `f - E x` in the next `m`.
         !! If `kode = 1` on entry, `res(i)` in {-1, 0, 1} for `i <= k` restricts the residual of the
         !! i-th equation to be <= 0, unrestricted or >= 0. Also used as scratch space
      real(wp), intent(out) :: error
         !! On exit the minimum sum of absolute values of the residuals `b - A x`

!     .. Workspace ..
      real(wp), allocatable :: cu(:, :)
         !! Cost coefficients, one row for each of the two signs of a variable
      integer, allocatable :: iu(:, :)
         !! Restriction flags, one row for each of the two signs of a variable
      integer, allocatable :: s(:)
         !! Row indices of the candidates of the ratio test
      character(len=256) :: alloc_msg
      integer :: alloc_stat
!     ..
!     .. Local Scalars ..
      real(wp) :: cuv, pivot, xmax, xmin, tmp1, zu, zv
      integer :: i, iq, ii, iineg, in, iout, iphase, j, js, kforce, kk
      integer :: klm, klm1, klm2, max_iter, n1, n2, nk
      logical :: to_phase2, at_optimum, pivot_found
!     ..
! CHECK THE ARGUMENTS BEFORE ANYTHING IS READ OR WRITTEN.
      if (k < 1 .or. l < 0 .or. m < 0 .or. n < 1 .or. &
          size(q, 1) < k + l + m + 2 .or. size(q, 2) < n + 2 .or. &
          size(x) < n .or. size(res) < k + l + m) then
         kode = L1_INVALID_INPUT
         iter = 0
         error = 0.0_wp
         return
      end if
! ALLOCATE THE WORKSPACE.
      allocate (cu(2, n + k + l + m), stat=alloc_stat, errmsg=alloc_msg)
      if (alloc_stat == 0) allocate (iu(2, n + k + l + m), stat=alloc_stat, errmsg=alloc_msg)
      if (alloc_stat == 0) allocate (s(k + l + m), stat=alloc_stat, errmsg=alloc_msg)
      if (alloc_stat /= 0) then
         kode = L1_ALLOC_FAILED
         iter = 0
         error = 0.0_wp
         return
      end if
!
! INITIALIZATION.
!
      max_iter = iter
      n1 = n + 1
      n2 = n + 2
      nk = n + k
      klm = k + l + m
      klm1 = klm + 1
      klm2 = klm + 2
      kforce = 1
      iter = 0
      js = 1
      iq = 0
      call set_up_labels(n, klm, q)
! SET UP PHASE 1 COSTS.
      call set_up_phase1_costs(n, k, l, m, q, cu, iu, iphase)
! NONNEGATIVITY RESTRICTIONS ON X (SIGN IN X(J)) AND ON THE RESIDUALS (SIGN IN RES(J)).
      if (kode /= 0) call apply_sign_restrictions(n, k, x, res, q, cu, iu, iphase)
!
! MAIN LOOP.  EACH PASS (RE)STARTS WITH THE PHASE 2 COSTS (ONLY WHEN TO_PHASE2 IS SET) AND
! THE MARGINAL COSTS, THEN PIVOTS UNTIL THE CURRENT PHASE IS OPTIMAL.
!
      to_phase2 = (iphase == 2)
      simplex: do
          if (to_phase2) then
! SET UP PHASE 2 COSTS.
              call set_up_phase2_costs(n, k, klm, iu, q, cu, iq, iphase)
          end if
! COMPUTE THE MARGINAL COSTS.
          call compute_marginal_costs(n, klm, js, cu, q)

          iterate: do
! DETERMINE THE VECTOR TO ENTER THE BASIS.
              xmax = 0.0_wp
              at_optimum = (js > n)
              if (.not. at_optimum) then
                  entering_candidates: do j = js,n
                      zu = q(klm1,j)
                      ii = int(q(klm2,j))
                      if (ii > 0) then
                          zv = -zu - cu(1,ii) - cu(2,ii)
                      else
                          ii = -ii
                          zv = zu
                          zu = -zu - cu(1,ii) - cu(2,ii)
                      end if
                      if (kforce == 1 .and. ii > n) cycle entering_candidates
                      candidate_u: block
                          if (iu(1,ii) == 1) exit candidate_u
                          if (zu <= xmax) exit candidate_u
                          xmax = zu
                          in = j
                      end block candidate_u
                      if (iu(2,ii) == 1) cycle entering_candidates
                      if (zv <= xmax) cycle entering_candidates
                      xmax = zv
                      in = j
                  end do entering_candidates
                  at_optimum = (xmax <= toler)
              end if
!
! TEST FOR OPTIMALITY.
!
              if (at_optimum) then
                  if (kforce == 0) then
                      if (iphase == 1) then
                          if (q(klm1,n1) <= toler) then
                              to_phase2 = .true.
                              cycle simplex
                          end if
                          kode = L1_INFEASIBLE
                      else
                          kode = L1_OPTIMAL
                      end if
                      exit simplex
                  end if
                  if (iphase == 1 .and. q(klm1,n1) <= toler) then
                      to_phase2 = .true.
                      cycle simplex
                  end if
                  kforce = 0
                  cycle iterate
              end if
              if (q(klm1,in) /= xmax) then
                  do i = 1,klm2
                      q(i,in) = -q(i,in)
                  end do
                  q(klm1,in) = xmax
              end if
!
! DETERMINE THE VECTOR TO LEAVE THE BASIS.
!
              pivot_found = .false.
              artificial_pivot: block
                  if (iphase == 1 .or. iq == 0) exit artificial_pivot
                  xmax = 0.0_wp
                  find_artificial_pivot: do i = 1,iq
                      tmp1 = abs(q(i,in))
                      if (tmp1 <= xmax) cycle find_artificial_pivot
                      xmax = tmp1
                      iout = i
                  end do find_artificial_pivot
                  if (xmax <= toler) exit artificial_pivot
                  call swap_rows(q, iq, iout, n2)
                  iout = iq
                  iq = iq - 1
                  pivot = q(iout,in)
                  pivot_found = .true.
              end block artificial_pivot

              if (.not. pivot_found) then
                  kk = 0
                  ratio_candidates: do i = 1,klm
                      tmp1 = q(i,in)
                      if (tmp1 <= toler) cycle ratio_candidates
                      kk = kk + 1
                      res(kk) = q(i,n1)/tmp1
                      s(kk) = i
                  end do ratio_candidates
! RATIO TEST: PICK THE SMALLEST RATIO, OR BYPASS THE VERTEX AND PICK AGAIN.
                  bypass: do
                      if (kk <= 0) then
                          kode = L1_ROUNDING_ERRORS
                          exit simplex
                      end if
                      xmin = res(1)
                      iout = s(1)
                      j = 1
                      if (kk /= 1) then
                          smallest_ratio: do i = 2,kk
                              if (res(i) >= xmin) cycle smallest_ratio
                              j = i
                              xmin = res(i)
                              iout = s(i)
                          end do smallest_ratio
                          res(j) = res(kk)
                          s(j) = s(kk)
                      end if
                      kk = kk - 1
                      pivot = q(iout,in)
                      ii = int(q(iout,n2))
                      if (iphase /= 1) then
                          if (ii < 0) then
                              iineg = -ii
                              if (iu(1,iineg) == 1) exit bypass
                          else
                              if (iu(2,ii) == 1) exit bypass
                          end if
                      end if
                      ii = abs(ii)
                      cuv = cu(1,ii) + cu(2,ii)
                      if (q(klm1,in)-pivot*cuv <= toler) exit bypass
! BYPASS INTERMEDIATE VERTICES.
                      do j = js,n1
                          tmp1 = q(iout,j)
                          q(klm1,j) = q(klm1,j) - tmp1*cuv
                          q(iout,j) = -tmp1
                      end do
                      q(iout,n2) = -q(iout,n2)
                  end do bypass
              end if
!
! GAUSS-JORDAN ELIMINATION.
!
              if (iter >= max_iter) then
                  kode = L1_MAX_ITERATIONS
                  exit simplex
              end if
              iter = iter + 1
              call pivot_tableau(n, klm, js, in, iout, pivot, q)
! THE LABEL OF THE VARIABLE THAT LEFT THE BASIS IS NOW IN Q(KLM2,IN).
              ii = int(abs(q(klm2,in)))
              if (iu(1,ii) == 0 .or. iu(2,ii) == 0) cycle iterate
              call swap_columns(q, in, js, klm2)
              js = js + 1
          end do iterate
      end do simplex
!
! PREPARE OUTPUT.
!
      call extract_solution(n, k, klm, q, x, res, error)

   end subroutine cl1

   subroutine set_up_labels(n, klm, q)
      !! Label the variables of the simplex tableau `q` and make every right-hand side non-negative.
      !!
      !! The unknown `x(j)` is labelled `j` in row `klm+2` of column `j`, and the residual variable
      !! of row `i` is labelled `n+i` in column `n+2` of row `i`. A row whose right-hand side
      !! `q(i,n+1)` is negative is multiplied by -1, including its label, so that a negative label
      !! marks a row (variable) that was flipped.
      integer, intent(in) :: n
         !! Number of columns of the matrices `A`, `C`, `E`
      integer, intent(in) :: klm
         !! Total number of rows `k+l+m` of the matrices `A`, `C`, `E`
      real(wp), contiguous, intent(inout) :: q(:, :)
         !! Simplex tableau with at least `klm+2` rows and `n+2` columns
      integer :: i, j

      do j = 1, n
         q(klm + 2, j) = real(j, wp)
      end do
      label_rows: do i = 1, klm
         q(i, n + 2) = real(n + i, wp)
         if (q(i, n + 1) >= 0.0_wp) cycle label_rows
         do j = 1, n + 2
            q(i, j) = -q(i, j)
         end do
      end do label_rows
   end subroutine set_up_labels

   ! The argument limit (6) is ignored for now; to be fixed together with the interface of cl1.
   ! allow(too-many-arguments)
   subroutine extract_solution(n, k, klm, q, x, res, error)
      !! Read the solution of the final simplex tableau `q`.
      !!
      !! The row `i` of `q` belongs to the variable whose label is stored in `q(i,n+2)`: a label
      !! `j <= n` is the unknown `x(j)`, a label `n+i` the residual `res(i)`. A negative label marks
      !! a variable whose sign was flipped. Variables that are not in the final basis are 0.
      !! The objective `error` is the sum of the residuals of the `k` equations `A x = b`.
      integer, intent(in) :: n
         !! Number of columns of the matrices `A`, `C`, `E`
      integer, intent(in) :: k
         !! Number of rows of the matrix `A`
      integer, intent(in) :: klm
         !! Total number of rows `k+l+m`
      real(wp), contiguous, intent(in) :: q(:, :)
         !! Final simplex tableau with at least `klm+2` rows and `n+2` columns
      real(wp), contiguous, intent(inout) :: x(:)
         !! Solution: `x(1:n)` is set, the other entries (`size(x) >= n`) are not touched
      real(wp), contiguous, intent(inout) :: res(:)
         !! Residuals: `res(1:klm)` is set, the other entries (`size(res) >= klm`) are not touched
      real(wp), intent(out) :: error
         !! Sum of the absolute values of the residuals `b - A x`
      real(dp) :: residual_sum
      real(wp) :: sign_of_variable
      integer :: i, label

      residual_sum = 0.0_dp
      x(1:n) = 0.0_wp
      res(1:klm) = 0.0_wp
      do i = 1, klm
         label = int(q(i, n + 2))
         if (label > 0) then
            sign_of_variable = 1.0_wp
         else
            label = -label
            sign_of_variable = -1.0_wp
         end if
         if (label <= n) then
            x(label) = sign_of_variable*q(i, n + 1)
         else
            res(label - n) = sign_of_variable*q(i, n + 1)
            if (label >= n + 1 .and. label <= n + k) residual_sum = residual_sum + real(q(i, n + 1), dp)
         end if
      end do
      error = real(residual_sum, wp)
   end subroutine extract_solution

   subroutine swap_rows(q, row1, row2, num_columns)
      !! Exchange the first `num_columns` entries of the rows `row1` and `row2` of `q`.
      !!
      !! Only the leading columns are exchanged: `q` may have more columns than the tableau uses.
      real(wp), contiguous, intent(inout) :: q(:, :)
         !! Simplex tableau
      integer, intent(in) :: row1
         !! First row
      integer, intent(in) :: row2
         !! Second row (may equal `row1`)
      integer, intent(in) :: num_columns
         !! Number of leading columns that are exchanged
      real(wp) :: tmp
      integer :: j

      do j = 1, num_columns
         tmp = q(row1, j)
         q(row1, j) = q(row2, j)
         q(row2, j) = tmp
      end do
   end subroutine swap_rows

   subroutine swap_columns(q, column1, column2, num_rows)
      !! Exchange the first `num_rows` entries of the columns `column1` and `column2` of `q`.
      !!
      !! Only the leading rows are exchanged: `q` may have more rows than the tableau uses.
      real(wp), contiguous, intent(inout) :: q(:, :)
         !! Simplex tableau
      integer, intent(in) :: column1
         !! First column
      integer, intent(in) :: column2
         !! Second column (may equal `column1`)
      integer, intent(in) :: num_rows
         !! Number of leading rows that are exchanged
      real(wp) :: tmp
      integer :: i

      do i = 1, num_rows
         tmp = q(i, column1)
         q(i, column1) = q(i, column2)
         q(i, column2) = tmp
      end do
   end subroutine swap_columns

   subroutine compute_marginal_costs(n, klm, first_column, costs, q)
      !! Compute the marginal costs and store them in row `klm+1` of the simplex tableau `q`.
      !!
      !! For each column `j` from `first_column` to `n+1` the entry `q(klm+1,j)` is the sum over the
      !! rows `i` of `q(i,j)` times the cost of the variable that is basic in row `i`. For the
      !! columns `first_column` to `n` the cost of the (non-basic) variable of the column itself is
      !! then subtracted. A negative label means that the sign of the variable was flipped, which
      !! selects the cost of the second sign (`costs(2,-label)`) instead of `costs(1,label)`.
      !! The sums are accumulated in the higher precision `dp`.
      integer, intent(in) :: n
         !! Number of columns of the matrices `A`, `C`, `E`
      integer, intent(in) :: klm
         !! Total number of rows `k+l+m`
      integer, intent(in) :: first_column
         !! First column of the tableau whose marginal costs are computed
      real(wp), contiguous, intent(in) :: costs(:, :)
         !! Costs of the variables, one row for each of the two signs of a variable
      real(wp), contiguous, intent(inout) :: q(:, :)
         !! Simplex tableau with at least `klm+2` rows and `n+2` columns; row `klm+1` is set
      real(dp) :: weighted_sum
      real(wp) :: variable_cost
      integer :: i, j, label

      do j = first_column, n + 1
         weighted_sum = 0.0_dp
         do i = 1, klm
            label = int(q(i, n + 2))
            if (label < 0) then
               variable_cost = costs(2, -label)
            else
               variable_cost = costs(1, label)
            end if
            weighted_sum = weighted_sum + real(q(i, j), dp)*real(variable_cost, dp)
         end do
         q(klm + 1, j) = real(weighted_sum, wp)
      end do
      do j = first_column, n
         label = int(q(klm + 2, j))
         if (label < 0) then
            variable_cost = costs(2, -label)
         else
            variable_cost = costs(1, label)
         end if
         q(klm + 1, j) = q(klm + 1, j) - variable_cost
      end do
   end subroutine compute_marginal_costs

   ! The argument limit (6) is ignored for now; to be fixed together with the interface of cl1.
   ! allow(too-many-arguments)
   subroutine set_up_phase1_costs(n, k, l, m, q, costs, restricted, phase)
      !! Set up the costs and the restrictions of the first phase of the simplex method.
      !!
      !! The variables of the linear program are numbered `1:n` for the unknowns `x`, `n+1:n+k` for
      !! the residuals of the `k` equations `A x = b`, `n+k+1:n+k+l` for the artificial variables of
      !! the `l` equality constraints and `n+k+l+1:n+k+l+m` for the slacks of the `m` inequality
      !! constraints. Each variable has two signs: `costs(1,j)` and `restricted(1,j)` belong to the
      !! negative sign, `costs(2,j)` and `restricted(2,j)` to the positive sign. A restricted sign
      !! (flag 1) is not allowed to enter the basis. The artificial variables of the equality
      !! constraints carry a cost of 1 for both signs and are restricted; the slacks of the
      !! inequality constraints carry a cost of 1 and are restricted for the positive sign only.
      integer, intent(in) :: n
         !! Number of columns of the matrices `A`, `C`, `E`
      integer, intent(in) :: k
         !! Number of rows of the matrix `A`
      integer, intent(in) :: l
         !! Number of rows of the matrix `C`
      integer, intent(in) :: m
         !! Number of rows of the matrix `E`
      real(wp), contiguous, intent(in) :: q(:, :)
         !! Simplex tableau with the labels set (see `set_up_labels`)
      real(wp), contiguous, intent(out) :: costs(:, :)
         !! Phase 1 costs of the `n+k+l+m` variables, one row for each sign
      integer, contiguous, intent(out) :: restricted(:, :)
         !! Restriction flags of the `n+k+l+m` variables, one row for each sign
      integer, intent(out) :: phase
         !! 1 if the first phase is needed (there are equality constraints, or an inequality
         !! constraint that the starting basis violates), 2 if the starting basis is feasible
      integer :: first_artificial, last_artificial, first_slack, last_slack

      first_artificial = n + k + 1
      last_artificial = n + k + l
      first_slack = last_artificial + 1
      last_slack = last_artificial + m

      phase = 2
      costs = 0.0_wp
      restricted = 0
      if (l /= 0) then
         costs(:, first_artificial:last_artificial) = 1.0_wp
         restricted(:, first_artificial:last_artificial) = 1
         phase = 1
      end if
      if (m /= 0) then
         costs(2, first_slack:last_slack) = 1.0_wp
         restricted(2, first_slack:last_slack) = 1
         if (any(q(k + l + 1:k + l + m, n + 2) < 0.0_wp)) phase = 1
      end if
   end subroutine set_up_phase1_costs

   ! The argument limit (6) is ignored for now; to be fixed together with the interface of cl1.
   ! allow(too-many-arguments)
   subroutine apply_sign_restrictions(n, k, x, res, q, costs, restricted, phase)
      !! Add the sign restrictions requested through `kode = 1` to the phase 1 costs.
      !!
      !! A restricted sign gets a cost of 1 and the restriction flag 1 (see `set_up_phase1_costs`).
      !! A restriction that the starting basis violates makes the first phase necessary.
      integer, intent(in) :: n
         !! Number of columns of the matrices `A`, `C`, `E`
      integer, intent(in) :: k
         !! Number of rows of the matrix `A`
      real(wp), contiguous, intent(in) :: x(:)
         !! Sign restriction of each unknown: -1 for `x(j) <= 0`, 0 for none, 1 for `x(j) >= 0`
      real(wp), contiguous, intent(in) :: res(:)
         !! Sign restriction of the residual of each of the `k` equations: -1 for `b - A x <= 0`,
         !! 0 for none, 1 for `b - A x >= 0`
      real(wp), contiguous, intent(in) :: q(:, :)
         !! Simplex tableau with the labels set (see `set_up_labels`)
      real(wp), contiguous, intent(inout) :: costs(:, :)
         !! Phase 1 costs, one row for each sign
      integer, contiguous, intent(inout) :: restricted(:, :)
         !! Restriction flags, one row for each sign
      integer, intent(inout) :: phase
         !! Set to 1 if a restriction is violated by the starting basis; otherwise unchanged
      integer :: j

      do j = 1, n
         if (x(j) < 0.0_wp) then
            costs(1, j) = 1.0_wp
            restricted(1, j) = 1
         else if (x(j) /= 0.0_wp) then
            costs(2, j) = 1.0_wp
            restricted(2, j) = 1
         end if
      end do
      do j = 1, k
         if (res(j) < 0.0_wp) then
            costs(1, n + j) = 1.0_wp
            restricted(1, n + j) = 1
            if (q(j, n + 2) > 0.0_wp) phase = 1
         else if (res(j) /= 0.0_wp) then
            costs(2, n + j) = 1.0_wp
            restricted(2, n + j) = 1
            if (q(j, n + 2) < 0.0_wp) phase = 1
         end if
      end do
   end subroutine apply_sign_restrictions

   ! The argument limit (6) is ignored for now; to be fixed together with the interface of cl1.
   ! allow(too-many-arguments)
   subroutine set_up_phase2_costs(n, k, klm, restricted, q, costs, num_restricted_rows, phase)
      !! Set up the costs of the second phase of the simplex method.
      !!
      !! In the second phase only the residuals of the `k` equations `A x = b` (the variables
      !! `n+1:n+k`) cost something. The variables that are basic at the end of the first phase and
      !! whose current sign is restricted (the artificial variables of the equality constraints, and
      !! the variables restricted through `kode = 1`) get a cost of 0, and their rows are moved to
      !! the top of the tableau, because they must leave the basis first.
      integer, intent(in) :: n
         !! Number of columns of the matrices `A`, `C`, `E`
      integer, intent(in) :: k
         !! Number of rows of the matrix `A`
      integer, intent(in) :: klm
         !! Total number of rows `k+l+m`
      integer, contiguous, intent(in) :: restricted(:, :)
         !! Restriction flags of the variables, one row for each sign
      real(wp), contiguous, intent(inout) :: q(:, :)
         !! Simplex tableau with at least `klm+2` rows and `n+2` columns: rows are exchanged
      real(wp), contiguous, intent(out) :: costs(:, :)
         !! Phase 2 costs of the variables, one row for each sign
      integer, intent(inout) :: num_restricted_rows
         !! Number of rows at the top of the tableau that hold a basic variable with a restricted
         !! sign: increased by one for each row that is moved
      integer, intent(out) :: phase
         !! Always 2
      integer :: i, label

      phase = 2
      costs = 0.0_wp
      costs(:, n + 1:n + k) = 1.0_wp
      move_restricted_rows: do i = 1, klm
         label = int(q(i, n + 2))
         if (label > 0) then
            if (restricted(1, label) == 0) cycle move_restricted_rows
            costs(1, label) = 0.0_wp
         else
            label = -label
            if (restricted(2, label) == 0) cycle move_restricted_rows
            costs(2, label) = 0.0_wp
         end if
         num_restricted_rows = num_restricted_rows + 1
         call swap_rows(q, num_restricted_rows, i, n + 2)
      end do move_restricted_rows
   end subroutine set_up_phase2_costs

   ! The argument limit (6) is ignored for now; to be fixed together with the interface of cl1.
   ! allow(too-many-arguments)
   subroutine pivot_tableau(n, klm, first_column, pivot_column, pivot_row, pivot_value, q)
      !! Perform one Gauss-Jordan pivot step on the simplex tableau `q`.
      !!
      !! The variable of the column `pivot_column` enters the basis in the row `pivot_row`: the
      !! pivot row is divided by the pivot, a multiple of it is added to every other row so that the
      !! pivot column becomes a unit vector (stored in the compact form of the revised tableau,
      !! with the pivot cell holding `1/pivot` and the rest of the column divided by `-pivot`),
      !! and the labels of the entering and the leaving variable are exchanged. Only the columns
      !! `first_column:n+1` and the rows `1:klm+1` are updated, the first columns having been
      !! retired earlier.
      !!
      !! After the call `q(klm+2,pivot_column)` holds the label of the variable that left the basis.
      !!
      !! NOTE: the original paper suggests replacing the column elimination below by a helper that
      !! adds a multiple of one column to another, for compilers that can pass a column of a 2-D
      !! array to a 1-D dummy argument; see `legacy/CALGO552.f`.
      integer, intent(in) :: n
         !! Number of columns of the matrices `A`, `C`, `E`
      integer, intent(in) :: klm
         !! Total number of rows `k+l+m`
      integer, intent(in) :: first_column
         !! First column of the tableau that is still active
      integer, intent(in) :: pivot_column
         !! Column of the entering variable
      integer, intent(in) :: pivot_row
         !! Row of the leaving variable
      real(wp), intent(in) :: pivot_value
         !! The pivot, `q(pivot_row, pivot_column)` on entry
      real(wp), contiguous, intent(inout) :: q(:, :)
         !! Simplex tableau with at least `klm+2` rows and `n+2` columns
      real(wp) :: multiplier, negative_pivot, label
      integer :: i, j

      do j = first_column, n + 1
         if (j /= pivot_column) q(pivot_row, j) = q(pivot_row, j)/pivot_value
      end do
      eliminate_columns: do j = first_column, n + 1
         if (j == pivot_column) cycle eliminate_columns
         multiplier = -q(pivot_row, j)
         do i = 1, klm + 1
            if (i /= pivot_row) q(i, j) = q(i, j) + multiplier*q(i, pivot_column)
         end do
      end do eliminate_columns
      negative_pivot = -pivot_value
      do i = 1, klm + 1
         if (i /= pivot_row) q(i, pivot_column) = q(i, pivot_column)/negative_pivot
      end do
      q(pivot_row, pivot_column) = 1.0_wp/pivot_value
      label = q(pivot_row, n + 2)
      q(pivot_row, n + 2) = q(klm + 2, pivot_column)
      q(klm + 2, pivot_column) = label
   end subroutine pivot_tableau

end module l1_calgo552
