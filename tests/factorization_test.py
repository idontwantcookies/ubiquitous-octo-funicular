import pytest


from src import factorization

@pytest.mark.parametrize("x,factors,phi", [
    [48, {2:4, 3:1}, 16],
    [11, {11:1}, 10]
])
def test_totient(x, factors, phi):
    assert factorization.totient(x, factors) == phi

@pytest.mark.parametrize('n', [12, 850903, 717967279050961])
def test_pollard_rho(n):
    x = factorization.pollard_rho_factor(n)
    y = n // x
    assert x * y == n

@pytest.mark.parametrize("n,f", [
    [12, {2: 2, 3: 1}],
    [717967279050961, {12657973: 1, 56720557: 1}],
    [100, {2: 2, 5: 2}]
])
def test_factors(n, f):
    assert factorization.pollard_rho_prime_power_decomposition(n) == f

@pytest.mark.parametrize("n,p,u,alpha", [
    [51, 3, 17, 1],
    [16, 2, 1, 4],
    [16, 3, 16, 0]
])
def test_factor_out(n, p, u, alpha):
    assert factorization.factor_out(n, p) == (u, alpha)

@pytest.mark.parametrize("n,primes,powers,u", [
    [2**3 * 3**2 * 11**3 * 13 * 17, [2, 3, 5, 7, 11], {-1:0, 2:3, 3:2, 5:0, 7:0, 11:3}, 13*17],
    [-22, [2, 3, 5], {-1:1, 2:1, 3:0, 5:0}, 11]
])
def test_factor_with_limited_primes(n, primes, powers, u):
    assert factorization.factor_with_limited_primes(n, primes) == (powers, u)
