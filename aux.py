from pprint import pp

from tp import factors, powmod, gcd, prime_miller_rabin


def totient(n:int):
    powers = factors(n)
    phi = 1
    for p, e in powers.items():
        phi *= p**(e - 1) * (p - 1)
    return phi

def subgroup(b:int, n:int):
    phi = totient(n)
    powers = []
    pi = 1
    for _ in range(1, phi + 1):
        pi *= b
        powers.append(pi % n)
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
