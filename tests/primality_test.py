import pytest

from src import primality
from big_numbers import primes


@pytest.mark.parametrize("n,P", [
    [10, [2, 3, 5, 7]],
    [30, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]],
    [5, [2, 3, 5]]
])
def test_erastosthenes_sieve(n, P):
    assert primality.eratosthenes_sieve(n) == P


@pytest.mark.parametrize("p,result", primes)
def test_primes_miller_rabin(p, result):
    result = bool(result)
    assert primality.prime_miller_rabin(p) == result

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
    assert primality.prime_miller_rabin(n) == result