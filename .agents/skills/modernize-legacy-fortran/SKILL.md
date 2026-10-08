---
name: modernize-legacy-fortran
description: Modernize legacy Fortran 77/90 source (fixed-form, GOTO, COMMON, EQUIVALENCE, implicit typing, assumed-size arrays, real*8/double precision) into clean free-form Fortran 2008/2013 modules following a strict style guide. Use when asked to refactor, port, upgrade, clean up or "modernize" old Fortran code (.f/.for/.f77/legacy .f90), e.g. classic netlib/ACM TOMS/CALGO routines.
---

# Modernize Legacy Fortran (F77/F90 -> F2013)

Act as an experienced modern-Fortran (2008/2013) programmer. Convert legacy code
into maintainable, free-form, module-based code **without changing numerical
behavior**. The conventions come from the community style guide at
<https://github.com/JorgeG94/fortran_programmer_llm> (condensed in
`references/style-guide.md`).

Supporting files (paths relative to this skill's directory):

- `references/legacy-patterns.md` - construct-by-construct "legacy -> modern" cookbook (GOTO, COMMON, EQUIVALENCE, DATA, ENTRY, alternate RETURN, statement functions, shared DO terminators, specific intrinsics, FORMAT, ...). **Read it before rewriting any code.**
- `references/fortitude.md` - how to use the **Fortitude** linter as the objective check for this style guide: commands, recommended `fortitude.toml`, rule-code <-> legacy-construct map, fix/suppression policy. **Read it before the first lint run.**
- `references/style-guide.md` - naming, required/forbidden/recommended practices, documentation style, common LLM mistakes. Read it before writing new code.
- `references/gpu.md` - OpenACC / CUDA Fortran notes. Only read if the user asks for GPU support or the code already contains such directives.

## Tooling: Fortitude is the gatekeeper

This skill relies on [Fortitude](https://github.com/PlasmaFAIR/fortitude) (`fortitude check`), a fast Fortran linter with
autofixes, to measure and enforce the style rules. Use it at every stage, not just at the end:

- Check availability first: `fortitude --version` (install with `uv tool install fortitude-lint@latest` or `pip install fortitude-lint`).
- Always run with preview rules on and a concise format, e.g. `fortitude check --preview --output-format=concise <files>`;
  use `fortitude explain <CODE>` to understand a rule before acting on it.
- Fortitude reads **free-form** Fortran. Convert fixed-form sources first (layer 1 below), then lint.
- Use `fortitude check --fix` (safe fixes) on *new* modernized files for mechanical issues; review the diff (`--diff`) and re-run the numerical comparison afterwards.
  Never use `--unsafe-fixes` without reviewing each change.
- A file is **not done** until `fortitude check --preview` reports no diagnostics for it (or each remaining one is justified by a documented `! allow(...)`/per-file ignore).
- Fortitude cannot verify numerical equivalence or naming conventions; the compile-and-compare step and the manual checklist are still required.

If Fortitude cannot be installed, say so, use the manual fallback greps in `references/fortitude.md`, and report that limitation.

## Core principle: preserve behavior first

Modernization is refactoring. Results must stay the same (to floating-point
tolerance) unless the user explicitly approves a change. Therefore:

1. **Never silently "fix" algorithms.** If you spot a bug, suspicious numerics
   or undefined behavior in the legacy code, keep the behavior, and report it
   (and optionally add a `! NOTE:` comment) instead.
2. **Keep the original file.** Write modernized code to a new file/location
   (or work under version control) so the original stays available for
   comparison. Do not delete legacy sources unless asked.
3. **Work in small, verifiable steps** and compile/test between steps.
4. **Do not change the public calling interface** unless asked. If the
   interface must change (e.g. grouping >6 arguments into a derived type),
   keep a thin backward-compatible wrapper or tell the user precisely what
   changed at call sites.

## Workflow

### 1. Survey
- Locate all legacy sources (`.f`, `.for`, `.f77`, `.f90`, `.F`), build files
  and any existing tests. Check for a build system (fpm.toml, CMake, Makefile)
  and which compiler is available (`gfortran`, `flang`, `ifx`, `nvfortran`).
- Inventory the constructs present, e.g. with `grep -inE`:
  `goto|common|equivalence|data |entry|external|implicit|\*8|double precision|dimension.*\(\*\)|assign |pause|if *\(.*\) *[0-9]+ *, *[0-9]+`.
- Build a call graph: which routines call which, what data is shared through
  COMMON / SAVE / DATA, what are the entry points used by callers (including
  non-Fortran callers, e.g. Python via f2py/ctypes - check names and
  argument order before changing them; `bind(C)` or f2py signature files may
  depend on them).

### 2. Establish a safety net (before touching code)
- Establish the lint baseline: `fortitude check --preview --output-format=concise` on the legacy sources that
  are already free-form (or after layer 1), and count diagnostics per rule code. Create/extend the project's `fortitude.toml`
  (template in `references/fortitude.md`) so the style-guide rules are enabled.
- Build the legacy code as-is and capture reference output for representative
  inputs (existing tests, or write a small driver). Prefer several cases,
  including edge cases (empty/degenerate sizes, ties, singular/infeasible
  problems, error return codes).
- Record the compiler and flags. Compile legacy code with warnings on:
  `gfortran -std=legacy -Wall -Wextra -fcheck=all -g` to surface latent bugs.
- Compare modernized results with a tolerance (relative `~1e-12` for double
  precision, or bit-for-bit where the algorithm is deterministic and
  operation order is unchanged). Report iteration counts / status codes too.

### 3. Transform in layers (compile and compare after each layer)
Apply in this order; each layer should be a reviewable step. After each layer: compile, run the reference comparison, and
re-run Fortitude to confirm the rule family targeted by that layer dropped to zero (rule codes are in brackets).

1. **Source form**: fixed-form -> free-form (`.f90`), continuation `&`,
   `!` comments, remove column-1 `C`/`*` comments, remove sequence numbers
   in columns 73-80, convert tabs, lower-case keywords. `fortitude check --fix` can then take care of keyword
   case/spacing, `::`, `[...]`, relational operators, named `end` statements, quotes, whitespace [S*, MOD011, MOD021].
2. **Declarations & types**: add `implicit none`; declare everything
   (derive types from the old implicit rules *i-n -> integer, otherwise real*,
   or from `IMPLICIT` statements); replace `real*8`, `double precision`,
   `real(8)` with `real(wp)`; suffix **all** real literals with `_wp`
   (`1.0d-6` -> `1.0e-6_wp`, and an un-suffixed `0.1` becomes `0.1_wp` only
   if the legacy code intended the single-precision value - see the
   cookbook); replace `integer*4` etc. with `int32` from `iso_fortran_env`.
   [C001, C021, C022, MOD001, MOD002, PORT011, PORT012, PORT021, OB031, OB061]
3. **Wrap into modules**: one module per file with project prefix, `private`
   by default, explicit `public` list, `use ..., only:`. Put procedures
   in `contains` so interfaces become explicit. Remove `external`
   declarations and old `real function` type prefixes on the function name.
   [C091, C092, C121, C122, C131, C132, S211, S212, S221]
4. **Arguments**: add `intent(in/out/inout)` to every dummy argument;
   replace assumed-size `a(*)` / explicit-size `a(n)` dummy arrays with
   assumed-shape `a(:)` / `a(:,:)` where it is safe (see the lower-bound and
   sequence-association caveats in the cookbook); drop dimension arguments
   (`n`, `lda`) that become redundant via `size()`, **unless** they carry
   semantic meaning (leading dimension of a sub-block, logical size < physical size). [C061, C071, C072]
5. **Control flow**: eliminate `GOTO`, arithmetic IF, computed GOTO, assigned
   GOTO, `PAUSE`, shared `DO` terminators, `DO` with `CONTINUE` labels and
   non-integer/real loop variables. Use `if/else`, `select case`, `do while`,
   `exit`/`cycle` with construct names, `block`, and early `return`.
   [OB041, OB051, OB081, OB091-OB094, C141-C143, C191, S251-S255]
6. **Data sharing**: replace `COMMON` / `EQUIVALENCE` / `BLOCK DATA` / `DATA`-
   initialized hidden state with module variables or (preferably) derived
   types passed explicitly. Replace `SAVE`d hidden state with explicit state types.
   [OB011, OB012, OB013, OB021, OB001, C081, C082, MOD051]
7. **Array & intrinsic idioms**: array syntax / `sum`, `maxval`, `minloc`,
   `merge`, `where`, `do concurrent` (only when iterations are provably
   independent), generic intrinsics instead of specific names
   (`dabs` -> `abs`, `dsqrt` -> `sqrt`, `dble` -> `real(..., wp)`, `float`, `idint`...).
   Be careful: whole-array rewrites may alter summation order and thus rounding.
8. **Structure & documentation**: named constants instead of magic numbers,
   `associate` for long expressions (watch the compiler caveats), `block` for
   scoping temporaries, FORD `!!` doc comments including units for physical
   quantities, `snake_case` names with verb+noun procedure names.
   [C031, C032, S901, S902; naming and docs are not checked by Fortitude - review manually]
9. **Error handling & I/O**: replace `STOP`/`PAUSE`/`print *` diagnostics in
   library routines with status codes (`integer, intent(out) :: stat`) or
   `error stop` with a message, and use the project's logging approach if one
   exists. Keep `write` unit numbers out of hard-coded literals (use
   `output_unit` / `error_unit` from `iso_fortran_env`, `newunit=` for files).
   [C032, C043, C181-C183, PORT001]
10. **Optional** (only if requested or clearly beneficial): type-bound
    procedures, `pure`/`elemental` marking, submodules, `allocatable`
    workspaces replacing work arrays passed by the caller, GPU directives.

### 4. Verify
- Run `fortitude check --preview --output-format=concise` on every modernized file: it must exit 0 (no diagnostics),
  or each remaining item must be a documented, justified suppression.
- Compile with strict flags, e.g.
  `gfortran -std=f2018 -Wall -Wextra -Wconversion -Wimplicit-interface -Wno-unused-dummy-argument -fcheck=all -fimplicit-none -g`.
  (`-std=f2018` is a superset of f2013; use `-std=f2008` if the project is pinned to it.)
  Resolve every warning you introduced; do not suppress warnings with flags to hide problems.
- Re-run the reference cases and compare against legacy output. Report the
  maximum absolute/relative difference.
- Run any existing project tests (`fpm test`, `ctest`, `pytest`, ...), and the
  wrapper/binding build if one exists.
- Review the non-linted checklist items by hand (naming, `pure`/side-effect-free functions, units/docs, single precision module).
- If you cannot compile or run (no compiler, missing inputs), say so
  explicitly. Never claim equivalence you did not verify.

### 5. Report
Summarize: files created/changed, Fortitude version + command + diagnostics before -> after (and any suppressions with reasons), constructs removed, interface changes,
anything intentionally left as-is and why, suspected latent bugs found in the
legacy code, and exactly what validation was run and its result.

## Final checklist (every modernized file must satisfy)

- [ ] `fortitude check --preview` is clean (exit 0) for the file, using the project's `fortitude.toml`.
- [ ] Free-form source, `.f90` (or `.F90` only if the preprocessor is needed); one module per file; no `.f`.
- [ ] Module has project prefix, `implicit none`, `private` default, explicit `public`s, and `use ..., only:` everywhere.
- [ ] Every dummy argument has `intent`; every procedure is in a module (explicit interface); no `external`.
- [ ] No `goto`, arithmetic `if`, computed/assigned `goto`, `common`, `equivalence`, `entry`, alternate `return`, `pause`, statement functions, assumed-size `(*)` arrays, `block data`, Hollerith.
- [ ] No `real*8`, `double precision`, `real(8)`; a single `wp` kind from one precision module, and **every** real literal carries `_wp`.
- [ ] Names: modules `prefix_name`, types `name_t`, variables/procedures `snake_case`, parameters `UPPER_CASE`, logical functions `is_/has_/can_`.
- [ ] Functions are side-effect free (`intent(in)` only); otherwise subroutine. `pure`/`elemental` where valid. No I/O inside `pure`.
- [ ] Public procedures take <= 6 arguments (group into derived types), except documented hot kernels.
- [ ] No unexplained magic numbers; units documented in comments for physical quantities.
- [ ] Nesting depth <= 3-4 (use `cycle`/early `return`); named loops when `exit`/`cycle` targets an outer loop.
- [ ] Prefer `allocatable` over `pointer`; allocatable strings `character(len=:), allocatable` for dynamic text.
- [ ] No hidden state: no implicit `save` (module variables / `data`-initialized locals) without explicit `save` or a state type.
- [ ] Behavior verified against the legacy implementation (or the lack of verification is clearly stated).

## Verification lessons (from modernizing CALGO 552)

- **Freeze oracles.** Keep the repaired original as a frozen oracle, and freeze each major stage you
  want to stay bit-identical to (for example the first double-precision free-form version before an
  interface change) in `legacy/`, excluded from linting. A refactor is then checked against an oracle of
  the same precision, bit-for-bit, on many generated inputs (compare every output with `np.array_equal`,
  not a tolerance).
- **Differential fuzzing needs targeted generators.** A general random generator missed three subtle
  slips; generators aimed at the odd branches (sign restrictions, exact-zero tolerance with small integer
  data, ties and duplicated rows) caught them. Check the corpus with mutation testing: plant a one-line
  bug (`<=` -> `<`, swapped branch, off-by-one size) and make sure the suite fails.
- **Never run memory-corrupting mutations in-process.** Undersized workspaces make the compiled code
  hang or abort. Run them only against the bounds-checked build (`-fcheck=all -O0`) in a child process,
  one mutation per run, with a time limit (macOS has no `timeout`; use the tool's timeout), and restore
  the source after each. Keep a backup copy of the file under test and verify it is restored.
- **Keep suspicious test failures separate from library failures.** Several failures here were wrong test
  data or wrong expected values (non-unique L1 optimum, uninitialised comparison arrays, array
  temporaries from non-`contiguous` helper dummies). Check the test before changing the library.
- **Validation that the C ABI cannot reach** (shapes derived from the dimension arguments) needs a small
  Fortran test program, built strictly and run from pytest.
- **Record gaps honestly** in the test README (what was mutation-tested, what was not).

## Things to avoid when modernizing

- Don't rename/reorder public routines or arguments without telling the user,
  especially when other languages call them.
- Don't convert loops to `do concurrent` or array expressions unless
  independence is certain; reductions and loop-carried dependencies are not allowed in `do concurrent`.
- Don't replace `pointer` with `allocatable` when aliasing is actually relied upon.
- Don't "optimize" or restructure the numerical algorithm as part of the
  modernization; propose that separately.
- Don't invent intrinsics/constants that don't exist: there is no built-in `pi`
  (define `PI = 4.0_wp*atan(1.0_wp)`), `random_number` is a subroutine,
  declarations must precede executable statements (use `block` for local scope),
  array constructors use `[ ]`, and never re-declare a variable (Fortran is case-insensitive).
- Don't introduce `print *` for diagnostics in library code.
