import pytest
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.linalg import GF2SolveUnavailable


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    """A test that needs the Rust solver for a system too large for the dense fallback is skipped, not failed,
    when `cargo` is not installed."""
    outcome = yield
    if outcome.excinfo is not None and issubclass(outcome.excinfo[0], GF2SolveUnavailable):
        outcome.force_exception(pytest.skip.Exception(f"needs the Rust solver (gf2solve): {outcome.excinfo[1]}"))


@pytest.fixture(scope="session")
def cx3():
    return SignedComplex(3)


@pytest.fixture(scope="session")
def cx4():
    return SignedComplex(4)


@pytest.fixture(scope="session")
def L2():
    return label_set(2)


@pytest.fixture(scope="session")
def L3():
    return label_set(3)


@pytest.fixture(scope="session")
def V32(cx3, L2):
    return violating_pairs(cx3, L2)


@pytest.fixture(scope="session")
def V33(cx3, L3):
    return violating_pairs(cx3, L3)
