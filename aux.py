from pprint import pp

from factorization import factors, gcd
from mod import powmod


def subgroup(b:int, n:int, phi:int):
    powers = []
    pi = 1
    for _ in range(1, phi + 1):
        pi = pi * b % n
        powers.append(pi)
        if pi == 1: break
    return powers

def is_generator(g, n, phi):
    if gcd(g, n) != 1: return False
    for p in factors(phi).keys():
        k = (phi) // p
        if powmod(g, k, n) == 1: return False
    return True
