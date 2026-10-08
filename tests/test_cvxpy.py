# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "cvxpy",
#     "numpy",
# ]
# ///

from typing import Literal

import cvxpy as cp
import numpy as np


def solve_l1_regression_as_lp_auto_reformulate(
    A: np.ndarray,
    b: np.ndarray,
    C: np.ndarray,
    d: np.ndarray,
    E: np.ndarray,
    f: np.ndarray,
    solver: Literal["HIGHS", "CBC", "SCS", "CLARABEL", "OSQP", "GLPK", "SCIP"],
):
    """Solves min_x ||b - Ax||_1 using linear programming via cvxpy"""

    num_vars = A.shape[1]
    K, N = A.shape
    L, N = C.shape
    M, N = E.shape
    print("Dimensions:")
    print(f"A: {K}x{N}")
    print(f"C: {L}x{N}")
    print(f"E: {M}x{N}")

    # ||b - Ax ||_1 = sum(| b_i - A[i, :] x_i|)

    # Define and solve the CVXPY problem.
    x = cp.Variable(num_vars)
    objective = cp.norm(b - A @ x, 1)
    constraints = [C @ x == d, E @ x <= f]
    prob = cp.Problem(cp.Minimize(objective), constraints=constraints)

    prob.solve(solver=solver, verbose=True)

    # Print the results
    print("Optimal Solution:")
    for i in range(num_vars):
        print(f"x[{i}] = {x.value[i]:.4g}")
    print(f"Objective Value = {prob.value}")
    print(prob.solver_stats)


if __name__ == "__main__":
    # A = np.array(
    #     [
    #         [2.0, 0.0, 1.0, 3.0, 1.0],
    #         [7.0, 4.0, 4.0, 15.0, 7.0],
    #         [9.0, 4.0, 7.0, 20.0, 6.0],
    #         [2.0, 2.0, 1.0, 5.0, 3.0],
    #         [9.0, 3.0, 2.0, 14.0, 10.0],
    #         [4.0, 5.0, 0.0, 9.0, 9.0],
    #         [4.0, 4.0, 9.0, 17.0, -1.0],
    #         [1.0, 6.0, 2.0, 9.0, 5.0],
    #     ]
    # )
    # b = np.array([7.0, 4.0, 7.0, 4.0, 0.0, 4.0, 9.0, 6.0])

    # C = np.array(
    #     [
    #         [0.0, 4.0, 5.0, 9.0, -1.0],
    #         [3.0, 2.0, 7.0, 12.0, -2.0],
    #         [3.0, 6.0, 12.0, 21.0, -3.0],
    #     ]
    # )

    # d = np.array([5.0, 1.0, 6.0])

    # E = np.array(
    #     [
    #         [0.0, 3.0, 6.0, 9.0, -3.0],
    #         [6.0, 2.0, 4.0, 12.0, 4.0],
    #     ]
    # )

    # f = np.array([5.0, 6.0])

    num_vars = 250
    np.random.seed(0)
    A = np.random.randn(1000, num_vars)
    b = np.random.randn(1000)
    C = np.random.randn(10, num_vars)
    d = np.random.randn(10)
    E = np.random.randn(200, num_vars)
    f = np.random.randn(200)

    solve_l1_regression_as_lp_auto_reformulate(A, f, C, d, E, f, solver="CARABEL")
