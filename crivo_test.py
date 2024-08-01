import pytest

import crivo


@pytest.mark.parametrize("n, x", [
    [5, 2],
    [4, 2],
    [10, 3],
    [1, 1]
])
def test_isqrt(n, x):
    assert crivo.isqrt(n) == x


@pytest.mark.parametrize("n,primes", [
    [10, [2, 3, 5, 7]],
    [30, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]],
    [5, [2, 3, 5]]
])
def test_erastosthenes_sieve(n, primes):
    assert crivo.eratosthenes_sieve(n) == primes


@pytest.mark.parametrize("n,s,t", [
    [2, 1, 1],
    [3, 0, 3],
    [16, 4, 1],
    [100, 2, 25]
])
def test_oddify(n, s, t):
    assert crivo.oddify(n) == (s, t)


@pytest.mark.parametrize("n,p,sq", [
    [4, 13, True],
    [2, 5, False],
    [3, 11, True]
])
def test_is_square(n, p, sq):
    assert crivo.is_square(n, p) == sq

@pytest.mark.parametrize("p", [13, 11, 7, 5, 3])
def test_find_non_square(p):
    x = crivo.find_non_square(p)
    for i in range(2, p):
        assert (i * i) % p != x

@pytest.mark.parametrize("a,p,d", [
    [2, 7, 3],
    [4, 7, 3],
    [5, 11, 2],
    [5, 11, 2],
])
def test_msqrt(a, p, d):
    sq = crivo.msqrt(a, p, d)
    assert (sq * sq) % p == a
