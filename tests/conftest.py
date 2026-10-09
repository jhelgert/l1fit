"""pytest configuration: makes the helper modules importable and lists the solvers under test."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from instances import load_instances
from solvers import available_solvers

_SOLVERS = available_solvers()
_INSTANCES = load_instances()


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    if "solver" in metafunc.fixturenames:
        metafunc.parametrize(
            "solver", list(_SOLVERS.values()), ids=list(_SOLVERS), scope="session"
        )
    if "inst" in metafunc.fixturenames:
        metafunc.parametrize("inst", _INSTANCES, ids=[i.name for i in _INSTANCES])


@pytest.fixture(scope="session")
def all_solvers():
    """All adapters by name (``src``: the Fortran sources, ``extension``: the installed package)."""
    return _SOLVERS
