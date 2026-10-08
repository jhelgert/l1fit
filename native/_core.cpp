// Thin nanobind binding of the Fortran L1 solver (ACM TOMS algorithm 552, see src/l1_calgo552.f90).
//
// The binding does no conversion at all. It validates the shapes of the arrays (the Fortran entry point
// derives them from the dimension arguments and cannot detect an inconsistent allocation), releases the
// GIL and calls the C entry point `l1_cl1` of src/l1_c_api.f90. The solver is reentrant (it has no module
// variables), so several threads may solve independent problems at the same time.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/tuple.h>

#include <cstdint>
#include <limits>
#include <string>
#include <tuple>

namespace nb = nanobind;

extern "C" {
// Fortran: bind(C, name="l1_cl1"). All arguments are passed by reference.
void l1_cl1(const int* k, const int* l, const int* m, const int* n, double* q, int* kode,
    const double* toler, int* iter, double* x, double* res, double* error);
}

namespace {

// Fortran (column-major) matrix and contiguous vectors of float64 in host memory. A non-const
// ndarray also rejects read-only arrays.
using Matrix = nb::ndarray<double, nb::ndim<2>, nb::f_contig, nb::device::cpu>;
using Vector = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw nb::value_error(message.c_str());
    }
}

std::tuple<int, int, double> cl1(int k, int l, int m, int n, Matrix q, int kode, double toler,
    int iter, Vector x, Vector res)
{
    require(k >= 0 && l >= 0 && m >= 0 && n >= 0, "the dimensions k, l, m, n must not be negative");
    // 64-bit arithmetic: k + l + m + 2 could overflow int
    const std::int64_t klm = static_cast<std::int64_t>(k) + l + m;
    require(klm + 2 <= std::numeric_limits<int>::max() && n + 2 <= std::numeric_limits<int>::max(),
        "the problem is too large");
    require(static_cast<std::int64_t>(q.shape(0)) == klm + 2 && static_cast<std::int64_t>(q.shape(1)) == static_cast<std::int64_t>(n) + 2,
        "q must have shape (k + l + m + 2, n + 2)");
    require(static_cast<std::int64_t>(x.shape(0)) == n, "x must have shape (n,)");
    require(static_cast<std::int64_t>(res.shape(0)) == klm, "res must have shape (k + l + m,)");

    double error = 0.0;
    {
        // The arguments are converted and checked: no Python object is touched from here on.
        nb::gil_scoped_release release;
        l1_cl1(&k, &l, &m, &n, q.data(), &kode, &toler, &iter, x.data(), res.data(), &error);
    }
    return { kode, iter, error };
}

} // namespace

NB_MODULE(_core, m)
{
    m.doc() = "Low-level binding of the Fortran L1 solver. Use l1fit.solve_l1 instead.";
    m.def("cl1", &cl1, nb::arg("k"), nb::arg("l"), nb::arg("m"), nb::arg("n"),
        nb::arg("q").noconvert(), nb::arg("kode"), nb::arg("toler"), nb::arg("iter"),
        nb::arg("x").noconvert(), nb::arg("res").noconvert(),
        "Solve min ||A x - b||_1 s.t. C x = d, E x <= f in place.\n\n"
        "`q` (Fortran-contiguous float64, shape (k+l+m+2, n+2)) holds [A b; C d; E f] and is destroyed,\n"
        "`x` (shape (n,)) and `res` (shape (k+l+m,)) are float64 vectors that receive the solution and\n"
        "the residuals (and, if kode == 1, hold sign restrictions on entry). The arrays are never\n"
        "converted or copied: a wrong dtype, layout or read-only array raises TypeError.\n"
        "Returns (status, iterations, objective); `iter` is the maximum number of iterations.");
}
