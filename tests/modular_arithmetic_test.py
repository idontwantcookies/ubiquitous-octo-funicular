import pytest
import sympy

from src import modular_arithmetic
from tests.big_numbers import modular_inverse, fastexp

@pytest.mark.parametrize("n,p,sq", [
    [4, 13, True],
    [2, 5, False],
    [3, 11, True]
])
def test_is_square(n, p, sq):
    assert modular_arithmetic.is_square(n, p) == sq

@pytest.mark.parametrize("p", [13, 11, 7, 5, 3])
def test_find_non_square(p):
    x = modular_arithmetic.find_non_square(p)
    for i in range(2, p):
        assert (i * i) % p != x

@pytest.mark.parametrize("a,p,d", [
    [2, 7, 3],
    [4, 7, 3],
    [5, 11, 2],
    [5, 11, 2],
])
def test_msqrt(a, p, d):
    sq = modular_arithmetic.msqrt(a, p, d)
    assert (sq * sq) % p == a

@pytest.mark.parametrize("b,e,n,x", fastexp)
def test_exp_binaria(b, e, n, x):
    assert modular_arithmetic.powmod(b, e, n) == x

@pytest.mark.parametrize("b,e,n,x", [
    [2, -1, 7, 4],
    [2, -2, 7, 2],
    [3, -1, 13, 9],
    [3, -2, 13, 3]
])
def test_powmod_negative_exponent(b, e, n, x):
    assert modular_arithmetic.powmod(b, e, n) == x


@pytest.mark.parametrize("a,n,inv", modular_inverse)
def test_invmod(a, n, inv):
    assert modular_arithmetic.invmod(a, n) == inv


@pytest.mark.parametrize("n", [11, 13, 17, 19, 101])
def test_find_generator(n):
    phi = n - 1
    f = sympy.factorint(phi)
    g = modular_arithmetic.find_generator(n, phi, f)
    for i in range(2, n - 1):
        assert modular_arithmetic.powmod(g, i, n) != 1

def test_congruence():
    r = [1, 2, 3]
    n = [5, 7, 11]
    x = modular_arithmetic.congruence_system(r, n)
    for i in range(3):
        assert x % n[i] == r[i]

@pytest.mark.parametrize("b,n,phi,G", [
    [2, 5, 4, [2, 4, 3, 1]],
    [4, 5, 4, [4, 1]]
])
def test_subgroup(b, n, phi, G):
    assert modular_arithmetic.subgroup(b, n, phi) == G

def test_congruence_system():
    a = [2, 3, 2]
    n = [3, 5, 7]
    x = modular_arithmetic.congruence_system(a, n)
    assert x == 23

@pytest.mark.parametrize("b,n,phi,F,expected", [
    [4, 5, 4, {2:2}, False],
    [2, 5, 4, {2:2}, True]
])
def test_is_generator(b, n, phi, F, expected):
    assert modular_arithmetic.is_generator(b, n, phi, F) == expected

@pytest.mark.parametrize("g,n,phi,F,o",[
    [2, 7, 6, {2:1, 3:1}, 3],
    [3, 7, 6, {2:1, 3:1}, 6],
    [3, 121, 110, {11:1, 2:1, 5:1}, 5],
    [5, 12, 4, {2:2}, 2],
    [13, 40, 12, {2:2, 3:1}, 4]
])
def test_order(g, n, phi, F, o):
    assert modular_arithmetic.order(g, n, phi, F) == o
