from collections import defaultdict

import pytest

from src import quadratic_sieve as qs

@pytest.mark.parametrize('ker,solutions',[
    [[[1,0,0], [0,1,0]], [[0,0,0], [0,1,0], [1,0,0], [1,1,0]]]
])
def test_kernel_solutions(ker, solutions):
    assert qs.kernel_solutions(ker) == solutions

@pytest.mark.parametrize('n,B', [
    [87463,43],
    [3000, 19],
    [100, 8],
    [10, 4],
    [3, 3]
])
def test_find_B(n, B):
    assert qs.find_B(n) == B

def test_join_powers():
    p1 = defaultdict(lambda: 0, {2: 3, 5: 1})
    p2 = defaultdict(lambda: 0, {3: 1, 5: 2})
    assert qs.join_powers(p1, p2) == {2: 3, 3: 1, 5: 3}

def test_euler_sieve_method():
    assert qs.euler_sieve_method(13, [2, 3, 5, 7, 11]) == [2, 3]

def test_setup():
    B, M, primes = qs.setup(100)
    assert (B, M) == (8, 9)
    assert primes == [-1, 2, 3, 7]

def test_build_matrix():
    primes = [2, 3, 5]
    S = {
        5: defaultdict(lambda: 0, {2: 3, 5: 1}),
        8: defaultdict(lambda: 0, {2: 1, 3: 2})
    }
    assert qs.build_matrix_of_powers(S, primes) == [
        [3, 1],
        [0, 2],
        [1, 0]
    ]

def test_compose():
    assert qs.compose({2: 3, 5: 1}) == 40
    assert qs.compose({-1: 1, 2: 3}) == -8

def test_compose_from_solution():
    S = {
        6: {2: 1, 3: 5},
        4: {5: 1, 11: 1},
        2: {5: 3, 7: 2, 11: 1}
    }
    a, b = qs.compose_from_solution(S, [0, 1, 1])
    assert a == 8
    assert b == 1925

def test_quadratic_sieve_aux():
    S = {}
    qs.quadratic_sieve_aux(17, 5, S, [-1, 2, 3, 5])
    assert S == {5: {-1:0, 2: 3, 3: 0, 5: 0}}

@pytest.mark.parametrize('n', [10, 50, 33, 100, 973, 1817, 2951, 8051, 87463, 10201030027])
def test_quadratic_sieve(n):
    d = qs.quadratic_sieve(n)
    assert d not in (1, n)
    assert n % d == 0
