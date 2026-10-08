# l1fit

A tiny Python package for blazingly fast curve fitting by an tailored primal
simplex algorithm and tailored dual simplex algorithms. In detail, this package
can solve the unconstrained/constrained optimization problems:

```math
\min_{x \in \mathbb{R}^n} \|Ax - b\|_1 \text{ s.t .} Cx = d, \;Dx \leq f, \\
```

with matrices $A \in \mathbb{R}^{m \times n}$, $C \in \mathbb{R}^{p \times n}$,
$D \in \mathbb{R}^{q \times n}$, $f \in \mathbb{R}^q$ and vectors
$b \in \mathbb{R}^m$, $d \in \mathbb{R}^p$ and $f \in \mathbb{R}^q$. Note
that the equality and inequality constraints are optional.

## TODO

- Create a mid-large set of benchmark problem instances, solve all of them via CVXPY
and sort them by size and solver status etc. such that we don't need to rely
on too simply and straightforward benchmark problems.
