from random import randint

from primality import prime_miller_rabin
from mod import invmod, powmod, gcd


def random_prime(bits:int=1024, max_attempts:int=1000):
    for _ in range(max_attempts):
        p = randint(0, 2**bits) | 1
        if prime_miller_rabin(p): return p
    else:
        raise ValueError("Failed to generate prime.")


def generate_keys(bits:int=1024) -> tuple[int, int, int]:
    p = random_prime(bits)
    q = random_prime(bits)
    n = p * q
    phi = (p - 1) * (q - 1)
    while True:
        e = randint(2, phi)
        if gcd(e, phi) == 1: break
    f = invmod(e, phi)
    return n, e, f


def encode(M: int, e: int, n: int) -> int:
    return powmod(M, e, n)

def decode(C: int, f: int, n: int) -> int:
    return powmod(C, f, n)
