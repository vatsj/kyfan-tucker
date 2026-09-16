import pytest
from kyfan import SignedComplex, label_set, violating_pairs


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
