# Fortitude: the linter for this workflow

[Fortitude](https://github.com/PlasmaFAIR/fortitude) (PlasmaFAIR) is a fast
Fortran linter (Rust, tree-sitter based, Ruff-inspired) with auto-fixes and
>100 rules. Docs: <https://fortitude.readthedocs.io>. It is the **primary
objective check** for the style-guide rules in this skill. It checks style and
correctness *statically*; it does **not** prove numerical equivalence - you still
need the compile + compare-against-legacy step.

## Setup

```bash
fortitude --version                      # check if installed (tested with 0.9.2)
uv tool install fortitude-lint@latest    # or: pip install fortitude-lint
```
If it cannot be installed (no network), say so, fall back to manual grep checks
(see the end of this file) and `gfortran` warnings, and state that in the report.

Rule codes below are from the stable docs at the time of writing; always
confirm against the installed version with `fortitude explain <CODE>` or
`fortitude explain --summary`, since preview rules and codes can change.

## Commands

```bash
fortitude check                                   # lint everything under cwd
fortitude check src/foo.f90                       # a file / dir / glob
fortitude check --output-format=concise src/      # one line per diagnostic
fortitude check --preview src/                    # include preview (unstable) rules, e.g. OB012, OB013, OB081, OB091-094, S271
fortitude check --select=C,OB,MOD,PORT,S --preview --output-format=concise src/foo.f90
fortitude check --fix src/foo.f90                 # apply safe automatic fixes
fortitude check --fix --preview src/foo.f90
fortitude check --diff src/foo.f90                # show what --fix would change without writing
fortitude explain C061 OB081                      # rationale + example for rules
fortitude explain obsolescent --summary           # whole category
fortitude config check                            # list config options
```

- Exit code: 0 = clean (or all fixed), 1 = violations remain, 2 = bad config/CLI/internal error. Use this as the
  "done" signal. Use `--statistics` (if available in the installed version)
  or count by code (`... --output-format=concise | awk '{print $2}' | sort | uniq -c | sort -rn`)
  to see which rules dominate and plan the next layer.
- **Fortitude parses free-form Fortran.** Fixed-form input yields a syntax error (E001) and all other
  violations are discarded for that file - this also happens for fixed-form code that merely has a `.f90` extension
  (e.g. `src/CALGO552.f90` in this repo). Legacy `.f`/`.for` files are not scanned by default (extension). Convert the *source form* first (layer 1 of the workflow), then lint.
- `--fix` modifies files in place and only applies *safe* fixes by default. Preview with `--diff`. Only run it on the new
  modernized file(s), never on the only copy of the legacy source. Do **not** use `--unsafe-fixes` blindly: review each
  change (e.g. C001 only adds `implicit none`, it does not declare the variables). After any `--fix`, **rebuild and
  re-run the reference comparison**; the literal/kind rewrites (C021, C022, MOD001, MOD002, PORT011/012/021) are the ones that can alter numerical results.
- To lint only what you changed, `--git-staged` or `--git-since <ref>` restrict checks to modified lines (fortitude >= 0.9.0).
- Suppress a justified diagnostic with an `allow` comment placed *before* the statement it covers (it applies to the whole
  next statement, e.g. a full module or procedure). Rule codes, names, or category names may be listed, comma separated:
  ```fortran
  ! allow(too-many-arguments)   ! hot kernel, kept flat for performance (see style guide)
  subroutine compute_flux_kernel(...)
  ```
  Unused/unknown/duplicate/disabled allow comments are themselves reported (FORT001-FORT005). Prefer fixing over suppressing;
  every suppression must be listed in the final report with its reason.

## Recommended configuration

Create `fortitude.toml` at the project root if none exists (or add the table to
`fpm.toml` under `[extra.fortitude.check]`, `pyproject.toml` under `[tool.fortitude.check]`).
The default rule set is small; this skill's style guide needs several off-by-default rules switched on:

```toml
[check]
preview = true
line-length = 132
target-std = "f2018"   # options: f95, f2003, f2008, f2018 (default), f2023; use "f2008" if the project is pinned to F2008
# Start from the default set plus the categories that implement the style guide
extend-select = [
  "C021", # no-real-suffix                  -> every real literal needs _wp
  "C022", # implicit-real-kind              -> real(wp), never bare `real`
  "C031", # magic-number-in-array-size      -> named constants
  "C032", # magic-io-unit                   -> output_unit / newunit
  "C043", # missing-action-specifier
  "C132", # default-public-accessibility    -> private by default
  "C142", # exit-or-cycle-in-unlabelled-loop
  "C183", # stat-without-message
  "MOD001", # double-precision              -> real(wp)
  "MOD002", # double-precision-literal      -> 1.0e-6_wp, not 1.0d-6
  "MOD031", # include-statement
  "PORT001", # non-portable-io-unit
  "S082", # multiple-statements-per-line
  "S102", "S103", "S104",   # whitespace layout
  "S201", # superfluous-implicit-none
  "S211", # multiple-modules (one module per file)
  "S212", # program-with-module
  "S221", # function-missing-result
  "S251", "S252", "S253", "S254", "S255", # useless return / superfluous else after return|cycle|exit|stop
  "S262", "S263",           # array-declaration style
  "S902", # too-many-arguments (style guide: <= 6)
  "S901", # too-complex (cyclomatic complexity; helps enforce "avoid deep nesting")
]

[check.complexity]
max-args = 6          # style guide: public procedures take <= 6 arguments
# max-complexity = 10 # tune to the code base; check `fortitude config check.complexity`

# Kernels that are explicitly allowed >6 args or high complexity should be
# listed per-file (or via an allow comment) instead of being globally ignored:
# [check.per-file-ignores]
# "src/abc_kernels.f90" = ["S902"]
```
Check that the installed version accepts every code (`fortitude explain <code>`) and option name
(`fortitude config check`, `fortitude config check.complexity`); drop or rename unknown ones rather than
working around errors. By default fortitude only scans modern extensions (`f90, F90, f95, ... pf`), so legacy `.f`/`.for`
files are not linted until converted to `.f90`.

If the project already has a fortitude/fpm/pyproject configuration, **extend**
it rather than replacing it, and tell the user what you added.

## Mapping: style-guide rule / legacy construct -> Fortitude rule

| Style-guide requirement or legacy construct | Rules | Auto-fix? |
|---|---|---|
| `implicit none` everywhere, no implicit typing | C001, C002, C003 | partially |
| `intent` on all dummy args | C061 | no (needs judgment) |
| No assumed-size `a(*)` | C071, C072 | no |
| No `external`; procedures in modules | C091, C092 | no |
| `use ..., only:`; intrinsic modules marked | C121, C122 | partly |
| `private` by default | C131, C132 | yes |
| No implicit save (initialization in declaration) | C081, C082, MOD051 | partly |
| Pointer components need default init | C101 | yes |
| Labeled `exit`/`cycle` in named loops; end labels | C141, C142, C143 | yes |
| Real literals carry a kind; no bare `real` | C021, C022 | yes (**verify numerics**) |
| No `double precision`, no `d` exponent | MOD001, MOD002 | yes (**verify numerics**) |
| No `real*8`, `real(8)`; no literal kind suffix `1.0_8` | PORT021, PORT011, PORT012 | partly |
| No magic numbers (array sizes, I/O units) | C031, C032, PORT001 | no |
| `[...]` array constructors | MOD011 | yes |
| `<`, `>=`, `==` instead of `.lt.`, `.ge.`, `.eq.` | MOD021 | yes |
| No `include` | MOD031, OB201, OB211 | partly |
| Single-line attribute declarations (no out-of-line `dimension`/`intent` statements) | MOD041 | yes |
| No statement functions | OB001 | - |
| No `common` / `equivalence` / `block data` | OB011, OB012, OB013 | no |
| No `entry` | OB021 | no |
| No specific intrinsic names (`dabs`, `dsqrt`, `amax1`, ...) | OB031 | yes |
| No computed GOTO / arithmetic IF / PAUSE | OB041, OB081, OB051 | no / no / no |
| Old character syntax `character*8`, `character*(*)` | OB061 | yes |
| No `forall` -> `do concurrent` | OB071 | partly |
| No labelled `do`, shared terminators, GOTO to `end do` | OB091, OB092, OB093, OB094 | no / partly |
| Unreachable code after `stop`/`return`/`exit` | C191 | no |
| Check `stat=` of `allocate`/`deallocate`/I/O | C181, C182, C183 | no |
| File extension `.f90`/`.F90` | S091 | no (rename) |
| Free-form hygiene: `::`, named `end` statements, double quotes, keyword case/space, whitespace | S071, S061, S241, S231-S233, S101, S002, S081 | yes |
| One module per file; no program+module in same file | S211, S212 | no |
| Functions with `result()` clause | S221 | yes |
| Early-return structure (avoid `else` after `return`/`cycle`/`exit`) | S252-S255, S251 | yes |
| Few arguments, low complexity | S902, S901 | no |
| Line length | S001 | no |

Not covered by Fortitude (verify by hand/compiler): numerical equivalence, naming conventions
(`_t` suffix, `UPPER_CASE` parameters, project prefix, `is_/has_` prefixes), pure/elemental
marking, units in comments, FORD docs, one precision module, functions without side effects, avoidance of
nested `associate`, magic numbers in executable code (only array sizes / I/O units are checked).

## How to use it in the loop

1. **Baseline** (after fixed -> free-form conversion): `fortitude check --preview --output-format=concise <file>`
   and record the count per rule code. This is the to-do list.
2. **After each layer**, re-run and confirm the relevant rule family went to zero
   (e.g. after "Control flow": OB041, OB081, OB091-094 -> 0; after "Declarations": C001, C021, C022, MOD001, MOD002, PORT011, PORT021 -> 0).
3. **Run `--fix` early** for mechanical rules (keyword case, `::`, `[...]`, relational operators, end
   labels, specific intrinsic names) to cut noise - then compile and diff results.
4. **Finish** with `fortitude check --preview` returning **no diagnostics** (exit 0) on every
   modernized file, then run the numerical comparison one last time since `--fix`
   may have touched code.
5. In the final report include: the exact fortitude command and version, number of diagnostics
   before -> after, and any allow-comments/ignores with reasons.

## Manual fallback checks (only if Fortitude is unavailable)

```bash
grep -inE '\bgo ?to\b|\bcommon\b|\bequivalence\b|\bentry\b|\bexternal\b|\bblock data\b|\bpause\b' src/*.f90
grep -inE 'double precision|real *\* *[0-9]|real *\( *[0-9]+ *\)|[0-9][dD][+-]?[0-9]' src/*.f90
grep -inE '\((\*|[a-z0-9_,: ]*,\*)\)' src/*.f90        # assumed-size
grep -inE '^ *use +[a-z_0-9]+ *$' src/*.f90            # use without only
```

## Quirks seen in practice

- **C061 false positive on a sibling's dummy name.** If a name is a *local* variable in one procedure of a
  module and a *dummy argument* of a sibling procedure (for example a local `cu` in the driver and a dummy
  `cu` in an extracted helper), Fortitude may report `C061 ... missing 'intent'` for the *local*
  declaration. Renaming the helper's dummy to a descriptive name (`costs`) fixes it, and is better style
  anyway. Extracting a block into a helper therefore needs distinct dummy names.
- **Unused `allow` comments are reported** (FORT002), so a suppression that stops being needed is flagged:
  useful to find every `allow(too-many-arguments)` that can be dropped once an interface is fixed.
