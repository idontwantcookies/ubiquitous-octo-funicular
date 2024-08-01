import pytest
import sympy

from base import isqrt, gcd_extended, poly, oddify
from modular_arithmetic import powmod, find_generator, congruence_system, order, invmod, is_square, find_non_square, msqrt, subgroup, is_generator
from primality import prime_miller_rabin, eratosthenes_sieve
from factorization import pollard_rho_factor, pollard_rho_prime_power_decomposition, totient
from discrete_log import baby_step_giant_step, pohlig_hellman, pohlig_hellman_prime_power_order


extgcd = []
with open("test/extended_gcd.txt") as file:
    for line in file:
        extgcd.append([int(x) for x in line.split()])

fastexp = []
with open("test/exp_binaria.txt") as file:
    for line in file:
        fastexp.append([int(x) for x in line.split()])

modular_inverse = []
with open("test/inverso_modular.txt") as file:
    for line in file:
        modular_inverse.append([int(x) for x in line.split()])

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
    assert isqrt(n) == r

@pytest.mark.parametrize("n,P", [
    [10, [2, 3, 5, 7]],
    [30, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]],
    [5, [2, 3, 5]]
])
def test_erastosthenes_sieve(n, P):
    assert eratosthenes_sieve(n) == P

@pytest.mark.parametrize("n,s,t", [
    [2, 1, 1],
    [3, 0, 3],
    [16, 4, 1],
    [100, 2, 25]
])
def test_oddify(n, s, t):
    assert oddify(n) == (s, t)

@pytest.mark.parametrize("n,p,sq", [
    [4, 13, True],
    [2, 5, False],
    [3, 11, True]
])
def test_is_square(n, p, sq):
    assert is_square(n, p) == sq

@pytest.mark.parametrize("p", [13, 11, 7, 5, 3])
def test_find_non_square(p):
    x = find_non_square(p)
    for i in range(2, p):
        assert (i * i) % p != x

@pytest.mark.parametrize("a,p,d", [
    [2, 7, 3],
    [4, 7, 3],
    [5, 11, 2],
    [5, 11, 2],
])
def test_msqrt(a, p, d):
    sq = msqrt(a, p, d)
    assert (sq * sq) % p == a


@pytest.mark.parametrize("a,b,x,y,d", extgcd)
def test_gcd_extended(a, b, x, y, d):
    assert gcd_extended(a, b) == (d, x, y)


@pytest.mark.parametrize("b,e,n,x", fastexp)
def test_exp_binaria(b, e, n, x):
    assert powmod(b, e, n) == x

@pytest.mark.parametrize("b,e,n,x", [
    [2, -1, 7, 4],
    [2, -2, 7, 2],
    [3, -1, 13, 9],
    [3, -2, 13, 3]
])
def test_powmod_negative_exponent(b, e, n, x):
    assert powmod(b, e, n) == x


@pytest.mark.parametrize("a,n,inv", modular_inverse)
def test_invmod(a, n, inv):
    assert invmod(a, n) == inv


@pytest.mark.parametrize("p,result", primes)
def test_primes_miller_rabin(p, result):
    result = bool(result)
    assert prime_miller_rabin(p) == result


@pytest.mark.parametrize("n,result", [
    [2, True],
    [3, True],
    [4, False],
    [10, False],
    [45, False],
    [101, True],
    [211, True],
    [21, False]
])
def test_miller_rabin_small(n, result):
    assert prime_miller_rabin(n) == result


@pytest.mark.parametrize("n", [11, 13, 17, 19, 101])
def test_find_generator(n):
    phi = n - 1
    f = sympy.factorint(phi)
    g = find_generator(n, phi, f)
    for i in range(2, n - 1):
        assert powmod(g, i, n) != 1


@pytest.mark.parametrize("g,x,h,p", bsgs)
def test_baby_step_giant_step(g, x, h, p):
    x = baby_step_giant_step(g, h, p, p - 1)
    assert x != 0 and x is not None
    assert g**x % p == h

@pytest.mark.parametrize("x,factors,phi", [
    [48, {2:4, 3:1}, 16],
    [11, {11:1}, 10]
])
def test_totient(x, factors, phi):
    assert totient(x, factors) == phi

def test_congruence():
    r = [1, 2, 3]
    n = [5, 7, 11]
    x = congruence_system(r, n)
    for i in range(3):
        assert x % n[i] == r[i]

@pytest.mark.parametrize("b,n,phi,G", [
    [2, 5, 4, [2, 4, 3, 1]],
    [4, 5, 4, [4, 1]]
])
def test_subgroup(b, n, phi, G):
    assert subgroup(b, n, phi) == G

def test_pohlig_hellman():
    n = 101
    f = {2: 2, 5: 2}
    g, h = 15, 100
    assert pohlig_hellman(g, h, n, f) == 50


def test_pohlig_hellman_prime_power():
    g, h, p, e, ord = 27, 40, 2, 3, 41
    assert pohlig_hellman_prime_power_order(g, h, p, e, ord) == 4


def test_poly():
    p = poly(1, 2, 3)
    assert p(2) == 11


@pytest.mark.parametrize('n', [12, 850903, 717967279050961])
def test_pollard_rho(n):
    x = pollard_rho_factor(n)
    y = n // x
    assert x * y == n


@pytest.mark.parametrize("n,f", [
    [12, {2: 2, 3: 1}],
    [717967279050961, {12657973: 1, 56720557: 1}],
    [100, {2: 2, 5: 2}]
])
def test_factors(n, f):
    assert pollard_rho_prime_power_decomposition(n) == f


def test_congruence_system():
    a = [2, 3, 2]
    n = [3, 5, 7]
    x = congruence_system(a, n)
    assert x == 23

@pytest.mark.parametrize("b,n,phi,F,expected", [
    [4, 5, 4, {2:2}, False],
    [2, 5, 4, {2:2}, True]
])
def test_is_generator(b, n, phi, F, expected):
    assert is_generator(b, n, phi, F) == expected

@pytest.mark.parametrize("g,n,phi,F,o",[
    [2, 7, 6, {2:1, 3:1}, 3],
    [3, 7, 6, {2:1, 3:1}, 6],
    [3, 121, 110, {11:1, 2:1, 5:1}, 5],
    [5, 12, 4, {2:2}, 2],
    [13, 40, 12, {2:2, 3:1}, 4]
])
def test_order(g, n, phi, F, o):
    assert order(g, n, phi, F) == o
