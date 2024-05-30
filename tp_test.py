import pytest
import sympy

import tp

extgcd = []
with open("test/extended_gcd.txt") as file:
    for line in file:
        extgcd.append([int(x) for x in line.split()])

fastexp = []
with open("test/exp_binaria.txt") as file:
    for line in file:
        fastexp.append([int(x) for x in line.split()])

invmod = []
with open("test/inverso_modular.txt") as file:
    for line in file:
        invmod.append([int(x) for x in line.split()])

bsgs = []
with open("test/bsgs.txt") as file:
    for line in file:
        bsgs.append([int(x) for x in line.split()])

primes = []
with open("test/primes.txt") as file:
    for line in file:
        primes.append([int(x) for x in line.split()])

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
    assert tp.isqrt(n) == r

@pytest.mark.parametrize("a,b,x,y,d", extgcd)
def test_gcd_extended(a, b, x, y, d):
    assert tp.gcd_extended(a, b) == (d, x, y)


@pytest.mark.parametrize("b,e,n,x", fastexp)
def test_exp_binaria(b, e, n, x):
    assert tp.powmod(b, e, n) == x


@pytest.mark.parametrize("a,n,inv", invmod)
def test_inverso_modular(a, n, inv):
    assert tp.invmod(a, n) == inv


@pytest.mark.parametrize("p,result", primes)
def test_primes_miller_rabin(p, result):
    result = bool(result)
    assert tp.prime_miller_rabin(p) == result

@pytest.mark.parametrize("n,result", [
    [2, True],
    [3, True],
    [4, False],
    [10, False],
    [101, True],
    [211, True],
    [21, False]
])
def test_miller_rabin_small(n, result):
    assert tp.prime_miller_rabin(n) == result


@pytest.mark.parametrize("n", [11, 13, 17, 19, 101])
def test_find_generator(n):
    factors = sympy.factorint(n - 1)
    g = tp.find_generator(n, factors)
    for i in range(2, n - 1):
        assert tp.powmod(g, i, n) != 1


@pytest.mark.parametrize("g,x,a,n", bsgs)
def test_baby_step_giant_step(g, x, a, n):
    x = tp.baby_step_giant_step(a, g, n)
    assert x != 0 and x is not None
    assert g**x % n == a


# @pytest.mark.parametrize("g,x,a,n", bsgs)
# def test_pohlig_hellman_base(g, x, a, n):
#     x = tp.pohlig_hellman_base(n, 1, g, a)
#     assert x != 0 and x is not None
#     assert g**x % n == a

def test_poly():
    p = tp.poly(1, 2, 3)
    assert p(2) == 11


@pytest.mark.parametrize('n', [12, 850903, 717967279050961])
def test_pollard_rho(n):
    x = tp.pollard_rho(n)
    y = n // x
    assert x * y == n


@pytest.mark.parametrize("n,factors", [
    [12, {2: 2, 3: 1}],
    [717967279050961, {12657973: 1, 56720557: 1}],
    [100, {2: 2, 5: 2}]
])
def test_factors(n, factors):
    assert tp.factors(n) == factors
