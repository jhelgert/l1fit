# Modern Fortran Style Guide (condensed)

Source: <https://github.com/JorgeG94/fortran_programmer_llm> (`Fortran_programmer.md`).
Target: modern Fortran (2008/2013+), emphasis on scientific computing.

## Naming conventions

| Item | Convention | Example |
|------|-----------|---------|
| Module file | `module_name.f90` lowercase with underscores; one module per file (submodules may be separate) | `swe_solver.f90` |
| Preprocessed file | `.F90` only when preprocessing is needed | `swe_config.F90` |
| Test file | `test_module_name.f90` | `test_swe_solver.f90` |
| Module | short 2-4 letter project prefix + snake_case | `swe_solver`, `swe_boundary` |
| Derived type | snake_case with `_t` suffix | `mesh_t`, `simulation_state_t` |
| Variable | descriptive snake_case | `num_cells`, `total_energy` |
| Procedure | verb + noun | `compute_flux`, `apply_boundary_conditions` |
| Logical function | `is_` / `has_` / `can_` prefix | `is_dry`, `has_converged` |
| Constant (parameter) | UPPER_CASE with underscores | `GRAVITY`, `MAX_ITERATIONS` |

- Bad: `ShallowWaterSolver`, `solver` (no prefix), `SimulationState`, `nCells`, `nc`, `E`, `flux(...)` (noun only), `maxIterations`.
- Single-letter names allowed only for loop indices (`i,j,k`), obvious local coordinates (`x,y,z`), and symbols matching published math notation.
- Choose one prefix and use it consistently so `use` statements are self-documenting and grep-able.
- For scientific code document physical units in declarations:
  ```fortran
  real(wp) :: depth        !! Water depth [m]
  real(wp), parameter :: GRAVITY = 9.81_wp  ! [m/s^2]
  ```

## Required practices

- **`use` with `only`**: `use iso_fortran_env, only: real64, int32`. Never bare `use netcdf`.
- **`implicit none`** in every module and program.
- **`intent`** for every dummy argument (`in`, `out`, `inout`).
- **Functions have no side effects**: all args `intent(in)`; use a subroutine to mutate arguments.
- **Private by default**: `private` at top of module, then explicit `public :: ...`.
- **<= 6 arguments** for public procedures; group related data into derived types (`mesh_t`, `state_t`, `config_t`, ...). Exception: stable hot compute kernels called millions of times may keep explicit scalar/array arguments to avoid derived-type indirection overhead.

## Forbidden practices

- `GOTO` -> structured control flow (`if`, `select case`, `exit`, `cycle`, `return`).
- Arithmetic `IF (x) 10,20,30` -> `if/else` or `select case`.
- `COMMON` blocks -> module variables or derived types.
- `EQUIVALENCE` -> proper conversion or `transfer()` if truly needed.
- Fixed-form source (`.f`, `.F`) -> free-form `.f90` / `.F90`.
- Assumed-size arrays `arr(*)` -> assumed-shape `arr(:)` with explicit interfaces (module procedures).
- `external` statements -> `use module, only: ...`.
- Implicit `save` (module variables with initializers, `data`-initialized locals) -> state in derived types passed as arguments; if unavoidable, write `save` explicitly.

## Recommended practices

- **Working precision**: one module defines `integer, parameter :: wp = real64`. Always `real(wp)`, literals `1.0e-10_wp`. Never `real(8)`, `real*8`, `double precision`.
- **Allocatable over pointer**. Use pointers only for aliasing existing data, linked structures, or polymorphic return.
- **`pure` / `elemental`** for side-effect-free procedures; elemental for scalar ops that should work on arrays.
- **`block`** constructs to limit variable scope and control deallocation timing of large workspaces.
- **`associate`** for long repetitive expressions. Caveats: support varies by compiler (aliasing/optimization problems in gfortran/nvfortran; flang most robust). Avoid *nested* `associate` constructs (known compiler problems). In performance-critical code fall back to explicit temporaries if problems appear.
- **No magic numbers**: named constants (`DRY_TOLERANCE`, `CONVERGENCE_TOL`); array sizes too.
- **Avoid deep nesting** (max 3-4 levels): invert conditions and use `cycle` / early `return`.
- **Labeled loops** for `cycle`/`exit` aimed at an outer loop (`outer: do ... exit outer`); optional but preferred over flag variables.
- **Allocatable strings**: `character(len=:), allocatable :: filename` rather than fixed `character(len=256)` when length is dynamic.
- **Memory management**: types with allocatable components are deallocated automatically; avoid needless `destroy` routines. Explicit cleanup is needed only for types holding pointers, GPU memory, and ordered resources (files, MPI communicators).
- **`do concurrent`** only for provably independent iterations; no `exit`, `cycle`, `return`, `goto` inside; compiler support varies; when in doubt use `do`.
- **Documentation**: FORD-compatible `!!` comments:
  ```fortran
  type :: simulation_t
     !! Main simulation container
     type(mesh_t) :: mesh
        !! Computational mesh
     real(wp) :: time = 0.0_wp
        !! Current simulation time [s]
  end type
  ```
- **Output**: use a logging facility rather than `print *`.

## Common LLM / AI mistakes to avoid

1. `pi` is not built in: `real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)`.
2. `random_number` is a subroutine: `call random_number(x)`.
3. No `print`/`write` inside `pure` procedures.
4. Declarations must come before executable statements (use `block` for local declarations mid-routine).
5. Don't declare a variable twice, including names differing only by case.
6. Array constructors: `[1, 2, 3]` or `(/ 1, 2, 3 /)`, never `(1, 2, 3)`.
7. Use named constants for array sizes in real code instead of literals.

## File structure template

```fortran
!! Brief module description (one line)
module abc_module_name
   !! Extended module documentation: purpose, usage, important notes.

   use iso_fortran_env, only: real64, int32
   use abc_other_module, only: needed_type_t, needed_function
   implicit none
   private

   ! Public API
   public :: my_type_t
   public :: initialize

   ! Constants
   real(real64), parameter :: SOME_CONSTANT = 1.0e-6_real64

   type :: my_type_t
      !! Type documentation
      integer :: n
         !! Number of elements
      real(real64), allocatable :: data(:)
         !! Data array
   contains
      procedure :: compute => my_type_compute
   end type my_type_t

contains

   subroutine initialize(self, n)
      !! Initialize the type with n elements
      type(my_type_t), intent(out) :: self
      integer, intent(in) :: n

      self%n = n
      allocate(self%data(n))
      self%data = 0.0_real64
   end subroutine initialize

   subroutine my_type_compute(self, input, output)
      !! Compute something
      class(my_type_t), intent(inout) :: self
      real(real64), intent(in) :: input
      real(real64), intent(out) :: output

      output = input * sum(self%data)
   end subroutine my_type_compute

end module abc_module_name
```

## Summary table

| Category | Do | Don't |
|----------|-----|-------|
| Types | `state_t`, `config_t` | `State`, `TConfig` |
| Variables | `num_cells`, `total_flux` | `nCells`, `tf` |
| Constants | `MAX_ITER`, `GRAVITY` | `maxIter`, `g` |
| Imports | `use mod, only: x, y` | `use mod` |
| Arrays | `arr(:)` | `arr(*)` |
| Memory | `allocatable` | `pointer` (unless needed) |
| Output | logging framework | `print *` |
| Control | `if` / `select case` / `do` | `goto` |

## Additional guidance from the Fortran Best Practices guide

Source: <https://fortran-lang.org/learn/best_practices/> (arrays, allocatable arrays, multidimensional
arrays, modules and programs, floating point, callbacks).

- **Arrays in procedures:** assumed-shape (`a(:)`, `a(:,:)`) is the default. The shape travels with the
  array and can be checked with `size()`/`shape()`. Explicit-shape arrays are unchecked (they can be passed
  with a wrong shape) and are meant for C interfaces and legacy code. Avoid passing whole slices (`f(r(:))`).
- **Contiguity:** assumed-shape dummies accept strided sections without a copy, so hot loops over them can
  be slower. Add `contiguous` to restore the unit-stride assumption of explicit-shape arrays; a
  non-contiguous actual argument is then copied in and out (visible with `-fcheck=array-temps`, part of
  `-fcheck=all`).
- **Scratch and work arrays:** use local `allocatable` arrays (heap) instead of caller-supplied workspace
  arguments (the "old LAPACK way") or automatic arrays that may overflow the stack. Allocate with
  `stat=`/`errmsg=`, one `allocate` per statement (Fortitude C181-C183), and report failure through a status.
- **Multidimensional arrays** are column-major: the inner loop should run over the first index.
  Row operations are inherently strided; do not restructure an algorithm for that.
- **Modules:** one module per file and the file name matches the module name; a library prefix; imports and
  `implicit none` at module scope; `private` with an explicit `public` list; documentation comment for the
  module, each procedure and the intent of each dummy argument. Avoid module variables (implicit `save`).
- **Kinds:** one central kind module; every literal carries a kind suffix. Integers convert exactly to
  `real(dp)` up to 2^53 (single precision only up to 2^24).
- **Callbacks / context:** prefer an abstract interface for procedure arguments and a derived type for
  context over `common` blocks, work arrays or `transfer`; `type(c_ptr)` is the equivalent of `void *`.
- **Names:** lowercase everywhere; short mathematical names (`k`, `n`, `q`, `x`) are acceptable for
  mathematical variables. This is looser than the "descriptive names" rule above: use it for the
  algorithm-internal variables and keep descriptive names for public arguments.
- **Layered interface for C / Python:** keep the Fortran API assumed-shape (shapes checked by the
  compiler) and put a thin `bind(C)` wrapper in its own module and file. The wrapper takes the
  dimensions, declares explicit-shape arrays, and forwards to the Fortran API (no copy). C has no shape
  information, so the wrapper cannot detect an inconsistent allocation: the caller (Python layer) must
  validate shapes. Do not use `error stop` for bad arguments in a library called from Python (it kills the
  interpreter); return a status code.
