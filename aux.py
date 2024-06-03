from pprint import pp

from tp import factors, powmod, gcd, totient


def subgroup(b:int, n:int):
    phi = totient(n)
    powers = []
    pi = 1
    for _ in range(1, phi + 1):
        pi = pi * b % n
        powers.append(pi)
        if pi == 1: break
    return powers

def is_generator(g, n):
    phi = totient(n)
    if gcd(g, n) != 1: return False
    for p in factors(phi).keys():
        k = (phi) // p
        if powmod(g, k, n) == 1: return False
    return True


for p in [3, 5, 7, 11, 13, 17, 23]:
    generators = []
    for x in range(2, p):
        if is_generator(x, p):
            generators.append(x)
    for e in range(2, 6):
        for g in generators:
            if not is_generator(g, p**e):
                print(f"{g} generates {p}, but does not generate {p}^{e}.")
                print(subgroup(g, p**e))
