# Legacy -> modern cookbook

Each section: what to look for, how to rewrite, and the numerical/semantic trap.
Fortitude rule codes that flag each construct are given in brackets (see `fortitude.md`).
`wp` is the working precision from a shared module, e.g.

```fortran
module abc_precision
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   integer, parameter, public :: wp = real64
end module abc_precision
```

---

## 1. Source form (fixed -> free) [S091, PORT031, E001 on unconverted input]

| Fixed form | Free form |
|-----------|-----------|
| `C` / `c` / `*` in column 1 | `!` |
| Continuation char in column 6 | trailing `&` on the previous line (optional leading `&`) |
| Statement label in columns 1-5 | label at start of line (only if still needed) |
| Columns 73-80 sequence numbers | delete (ignored in fixed form; in free form they would become code) |
| Blanks insignificant (`DO 10 I = 1,N`, `GO TO`) | blanks are significant: normalize `GOTO`/`GO TO`, `ENDIF` -> `end if` |
| Tabs | spaces |
| `INCLUDE 'file.inc'` [MOD031] | module (COMMON/PARAMETER includes become modules or types) |

Also: `.f` -> `.f90`; line length <= 132 (or the project's configured limit); 3-space indentation; lower-case keywords; double quotes for strings [S241]; named `end` statements [S061].

## 2. Implicit typing [C001, C002, C003]

Add `implicit none`, then declare every identifier. Find undeclared names by compiling the
legacy source with `gfortran -fimplicit-none -fsyntax-only`. Default rule: `a-h, o-z` -> `real`,
`i-n` -> `integer`, unless an `IMPLICIT` statement says otherwise. Beware a default-`real`
(single precision!) variable in a `double precision` program: it may be an accident or intentional - check how it is used,
and flag it in the report if converting it to `real(wp)` changes results.

## 3. Kinds and literals [MOD001, MOD002, C021, C022, PORT011, PORT012, PORT021, OB031, OB061]

| Legacy | Modern |
|--------|--------|
| `double precision x`, `real*8 x`, `real(8) x`, `real(kind=8) x` | `real(wp) :: x` |
| plain `real x` in a "double" code | decide: `real(wp)` if working precision was intended; `real(real32)` if single is intentional |
| `integer*4` / `integer*8` | `integer(int32)` / `integer(int64)` |
| `complex*16` | `complex(wp)` |
| `character*8 s`, `character*(*) s` | `character(len=8) :: s`, `character(len=*), intent(in) :: s` |
| `1.0D0`, `1.D-6`, `2.5d3` | `1.0_wp`, `1.0e-6_wp`, `2.5e3_wp` |
| `1.`, `.5` | `1.0_wp`, `0.5_wp` [S291] |
| `0.1` (no suffix) in double-precision expressions | `0.1_wp` - **trap**: legacy value is the single-precision 0.1 widened (`0.100000001490116`). Changes results at ~1e-8 relative level. Usually `0.1_wp` was intended; report it. Exactly representable values (`0.5`, `2.0`) are safe. |
| `DBLE(x)`, `SNGL(x)`, `FLOAT(i)`, `DFLOAT(i)` | `real(x, wp)`, `real(x)`, `real(i, wp)` |
| `DABS DSQRT DEXP DLOG DSIGN DMOD DMAX1 DMIN1 AMAX1 AMIN1 MAX0 MIN0 IDINT IABS` | `abs sqrt exp log sign mod max min max min int abs` (generic) - **trap**: `max(a, b)` with mixed kinds needs matching kinds |
| `D1MACH`, `R1MACH`, `I1MACH`, hand-rolled epsilon loops | `epsilon(1.0_wp)`, `tiny(1.0_wp)`, `huge(1.0_wp)`, `huge(1)`, `radix`, `digits` |
| `.LT. .LE. .EQ. .NE. .GT. .GE.` [MOD021] | `< <= == /= > >=` |
| `(/ 1, 2 /)` [MOD011] | `[1, 2]` |

Define `PI = 4.0_wp*atan(1.0_wp)`; there is no built-in `pi`.
Only introduce `ZERO`/`ONE` parameters if it improves readability; do not obscure the math.

## 4. GOTO and friends [OB041, OB081, OB051, OB094, C191]

### Forward GOTO / error exit
```fortran
      IF (IERR .NE. 0) GO TO 900
      ...
  900 CONTINUE
      RETURN
```
->
```fortran
if (ierr /= 0) then
   stat = ierr
   return
end if
```
For cleanup shared by several exits use a named `block` and `exit`:
```fortran
main: block
   if (cond_a) exit main
   ...
   if (cond_b) exit main
   success = .true.
end block main
if (.not. success) call report_failure()
```

### Backward GOTO (hand-made loop)
```fortran
   10 CONTINUE
      ... work ...
      IF (.NOT. DONE) GO TO 10
```
-> `do ... if (done) exit ... end do` (bottom-tested) or `do while (cond)` (top-tested). Check where the test is
evaluated relative to the work, and what values the variables hold at loop exit.

### Jump into the middle of a loop / spaghetti
Draw the control-flow graph first. Options: a logical flag; an integer `state` with `select case` inside a
`do` loop (`state = NEXT; cycle`); or split chunks into internal procedures. Preserve evaluation order exactly.
"Restart from step N" semantics (typical in simplex/pivot codes) -> outer named `do` with `cycle outer`.

### `GO TO` an `end do` / `continue` label inside a loop [OB094]
That is a `cycle`:
```fortran
      DO 20 I = 1, N
         IF (X(I) .LE. 0) GO TO 20
         ...
   20 CONTINUE
```
->
```fortran
do i = 1, n
   if (x(i) <= 0.0_wp) cycle
   ...
end do
```

### Arithmetic IF [OB081]
`IF (X) 10, 20, 30` ->
```fortran
if (x < 0) then
   ! 10
else if (x == 0) then
   ! 20
else
   ! 30
end if
```
Two-label forms (`IF (X) 10,20,20`) simplify to a plain `if/else`. NaN falls to the last branch in most
compilers; the `else` form behaves the same - mention if the legacy code depends on it.
Fall-through labels (a label that is the next statement) mean that branch is empty.

### Computed GOTO [OB041]
`GO TO (10,20,30), K` -> `select case (k)`; `case (1) ... case (2) ... case (3) ...`. When K is out of
range legacy code falls through to the next statement, so `case default` should be a no-op (not an error)
unless you verified otherwise [C011 wants a default case].

### Assigned GOTO / `ASSIGN`
Replace with an integer state variable + `select case`, or a procedure pointer / internal procedure.

### `PAUSE` [OB051], `STOP`
Delete `PAUSE`. In library routines replace `STOP` with a status argument; in programs use `stop` / `error stop "message"`.

### Alternate RETURN (`SUBROUTINE F(A, *, *)`, `RETURN 1`, `CALL F(A, *10, *20)`)
Add `integer, intent(out) :: stat` and `select case (stat)` at call sites.

## 5. DO loops [OB091, OB092, OB093]

| Legacy | Modern |
|--------|--------|
| `DO 10 I=1,N ... 10 CONTINUE` | `do i = 1, n ... end do` |
| Shared terminator: `DO 10 I=..` / `DO 10 J=..` / `10 CONTINUE` | two separate `do`/`end do` (inner first) |
| Terminating on an executable statement (`DO 10 I=1,N` / `10 A(I)=0`) | put the statement inside `do ... end do` |
| Real loop variable `DO 10 X=0.0,1.0,0.1` | integer-count loop; compute `x = x0 + k*dx`; trip count is `max(int((stop-start+step)/step),0)` - **verify the count** |
| Loop variable modified in the body | restructure into `do while` / explicit `exit` |
| Loop index used after the loop | value is `last + step` after normal completion; keep or assign explicitly |
| Nested loops that `exit`/`cycle` an outer loop | name the loops (`outer: do`) [C141, C142] |

Loop order: Fortran is column-major - innermost loop should run over the first index. Interchange only if
results (summation order) are allowed to change.

## 6. COMMON blocks [OB011] and BLOCK DATA [OB013]

1. Collect every `COMMON /name/` occurrence across all files. Compare layouts: legacy code often redeclares
   the same block with *different* names/shapes in different routines (storage association) - hidden aliasing to resolve deliberately.
2. Consistent layouts -> derived type `xxx_t` passed explicitly (`intent(in)` when read-only). Truly global singleton state ->
   module variable with explicit `save`, as a last resort.
3. Inconsistent layouts (storage reinterpretation) -> one owner array plus `associate` views / pointer rank remapping /
   explicit copies; document each reinterpretation.
4. Blank common (`COMMON // WORK(1000)`) is usually scratch -> `allocatable` local/workspace arrays.
5. `BLOCK DATA` initializing COMMON -> default component initialization in the type or a `parameter`/init routine.
6. Report each routine that gained a derived-type argument; keep a compatibility wrapper if external callers exist.

```fortran
module abc_grid
   use abc_precision, only: wp
   implicit none
   private
   public :: grid_t

   type :: grid_t
      integer  :: nx = 0
      integer  :: ny = 0
      real(wp) :: dx = 0.0_wp
      real(wp) :: dy = 0.0_wp
   end type grid_t
end module abc_grid
```

## 7. EQUIVALENCE [OB012]

- Memory saving / scratch reuse -> separate `allocatable` arrays, or `block` scoping.
- Aliasing a scalar and an array element -> use the element directly, or `associate`.
- Type punning (real <-> integer bits) -> `transfer()`; document it (processor-dependent as before).
- 2-D array viewed as 1-D -> `reshape()` copy, or rank-remapping pointer `p(1:n*m) => a` (`a` needs `target`).

## 8. DATA, SAVE and hidden state [C081, C082, MOD051]

- Variables with `DATA` or an initializer in the declaration (`integer :: k = 0`) are **implicitly saved**. Decide:
  (a) constant -> `parameter`; (b) true persistent state -> state type passed in, or at minimum an explicit `save`;
  (c) meant to reset each call -> declare plain and assign in the executable part (this changes behavior if the legacy code really relied on the saved value - check).
- Large `DATA` tables -> `parameter` arrays with constructors: `real(wp), parameter :: COEF(3) = [1.0_wp, 2.0_wp, 3.0_wp]`; 2-D via `reshape` (both fill column-major).
- `DATA` with implied-DO -> `[(expr, i = 1, n)]`.
- Old compilers kept all locals static; code that reads a local before setting it may "work" only by accident. Detect with
  `-finit-real=snan -ffpe-trap=invalid -fcheck=all`, report findings.
- `save` state makes the routine non-reentrant / not thread-safe; mention it.
- Pointer components need a default initializer `=> null()` [C101].

## 9. ENTRY [OB021], statement functions [OB001], EXTERNAL [C091, C092]

- `ENTRY` -> separate procedures sharing a private helper.
- Statement function `F(X) = X*X + A` -> internal (`contains`) `pure` function (host association gives access to `A`) or pass `a` explicitly.
- `EXTERNAL f` / `REAL f` passed as an argument -> abstract interface:
```fortran
abstract interface
   pure function objective_i(x) result(fx)
      import :: wp
      real(wp), intent(in) :: x(:)
      real(wp) :: fx
   end function objective_i
end interface
...
procedure(objective_i) :: fun     ! dummy procedure
```
- Every procedure goes in a module (`contains`) so interfaces are explicit; function declarations use `result()` [S221].

## 10. Arrays and argument passing [C071, C072, C061]

### Assumed-size and explicit-shape dummies
```fortran
      SUBROUTINE SOLVE(A, LDA, N, B)
      DOUBLE PRECISION A(LDA,*), B(N)
```
->
```fortran
subroutine solve(a, b)
   real(wp), intent(inout) :: a(:,:)
   real(wp), intent(inout) :: b(:)
```
Check each of these before converting:

- **Lower bounds**: assumed-shape dummies start at 1 regardless of the actual argument. For `A(0:N)` write `a(0:)`.
- **Sequence association**: legacy code passes `A(1,J)` (an element) to a `X(*)` dummy to mean "column J from here", passes a 2-D array to a 1-D dummy,
  or slices one work vector into several arrays by offset (`CALL F(W(I1), W(I2))`). Assumed-shape forbids this. Rewrite call sites with
  slices (`a(:,j)`, `work(i1:i1+n-1)`) or keep an explicit-shape dummy `a(n)` where reinterpretation is intended.
- **Leading dimension != extent** (`LDA > M`): pass the sub-array `a(1:m,1:n)`; `lda` disappears.
- **Logical size < physical size**: keep an explicit `n` argument (don't use `size()`) when the routine only works on the first `n` elements.
- **Non-contiguous slices** (`a(i,:)`) may trigger copy-in/out temporaries on hot paths; reorder loops or pass contiguous columns.
- `character(len=*)` dummies must have `intent(in)` unless they are really modified [C072].

### Intent
Infer `intent` from usage: only read -> `in`; assigned before read -> `out`; both -> `inout`. Trap: an argument assigned in only
some branches is not `out` (it keeps its old value otherwise) - use `inout`. `intent(out)` also resets default-initialized components and deallocates allocatables.

### Work arrays
Legacy `SUBROUTINE F(A, N, WORK, LWORK)` -> remove `work`/`lwork` and use local `allocatable` arrays (`allocate(work(n))`) unless
the caller needs to control memory; add `stat=`/`errmsg=` on large allocations [C181, C182].

### Too many arguments [S902]
When a public routine has > 6 arguments, bundle related ones into derived types (problem data, options, workspace/state, results).
Keep hot numerical kernels flat (explicit arrays) and suppress S902 for them with a reasoned `! allow(too-many-arguments)`.

## 11. Intrinsics, array idioms and `do concurrent` [OB071]

- Replace hand-written loops with `sum`, `maxval`, `minloc`, `dot_product`, `matmul`, `merge`, `where`, `count`, `any/all`, array sections - **but** a changed summation order changes rounding; keep explicit loops in code where bitwise reproducibility matters.
- `FORALL` -> `do concurrent`. `do concurrent` only when iterations are independent (no reductions or loop-carried dependencies; no `exit`/`cycle`/`return`/`goto` inside).
- Initialization `DO I=1,N; A(I)=0.0; END DO` -> `a = 0.0_wp`.
- Swap loops -> `call swap` via temporary or array assignment `tmp = a; a = b; b = tmp`.
- `MIN`/`MAX` with an index search -> `minloc`/`maxloc` (note they return the *first* occurrence; check tie-breaking matches the legacy loop).
- `ISIGN`/`SIGN` with signed zero and `abs`: keep the same intrinsic - semantic differences exist only for -0.0 and NaN.

## 12. I/O, STOP, errors [C032, PORT001, C043, C181, C183]

- Hard-coded units (`WRITE(6,...)`, `WRITE(5,...)`, `READ(5,...)`) -> `output_unit`, `error_unit`, `input_unit` from `iso_fortran_env`. Other
  fixed unit numbers -> `open(newunit=u, ..., action=..., iostat=ios, iomsg=msg)`; always check `iostat`.
- `FORMAT` statements -> inline format strings or character parameters (`character(len=*), parameter :: FMT_RESULT = "(a, i0, 2x, es12.4)"`); `Hollerith` constants (`4HABCD`) -> string literals.
- Library routines should return a status instead of printing/stopping. For diagnostics in programs use the project's logging approach, else `write(error_unit, "(a)")`. Avoid `print *`.
- Never put I/O in `pure` procedures.
- `ERR=`/`END=` labels on I/O -> `iostat=` + `if (ios /= 0)` / `is_iostat_end(ios)`.

## 13. Magic numbers and constants [C031, C032]

Replace unexplained literals (tolerances, iteration limits, array sizes, flag codes, return codes) with named `UPPER_CASE` parameters.
Algorithmic constants that come from a paper (e.g. `0.5`, `2.0`) can stay inline when they are part of an obvious formula.
Return/status codes: define them as parameters (`integer, parameter :: STATUS_OPTIMAL = 0, STATUS_UNBOUNDED = 2`) and document them.

## 14. Naming and documentation

- Rename terse legacy names (`NCT`, `IA`, `XX`) to descriptive snake_case **only inside the routine** where the meaning is certain; otherwise
  keep the original name and add a `!!` doc comment with the meaning (and units) [see style-guide.md]. Keep a mapping table old -> new in the report for renamed *public* items.
- Add a FORD `!!` header to each module and procedure: purpose, argument meanings, algorithm reference (paper, ACM TOMS number, authors, license).
- Preserve credits/license/citation text from the original file header.
