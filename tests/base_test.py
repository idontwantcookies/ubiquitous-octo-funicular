import pytest

import base
from tests.big_numbers import extgcd


@pytest.mark.parametrize("n,r", [
    [1, 1],
    [2, 1],
    [3, 1],
    [4, 2],
    [5, 2],
    [6, 2],
    [7, 2],
    [8, 2],
    [9, 3],
    [101, 10]
])
def test_isqrt(n, r):
    assert base.isqrt(n) == r

@pytest.mark.parametrize("n,s,t", [
    [2, 1, 1],
    [3, 0, 3],
    [16, 4, 1],
    [100, 2, 25]
])
def test_oddify(n, s, t):
    assert base.oddify(n) == (s, t)

@pytest.mark.parametrize("a,b,x,y,d", extgcd)
def test_gcd_extended(a, b, x, y, d):
    assert base.gcd_extended(a, b) == (d, x, y)

def test_poly():
    p = base.poly(1, 2, 3)
    assert p(2) == 11
