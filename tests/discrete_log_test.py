import pytest

from src import discrete_log
from tests.big_numbers import bsgs


@pytest.mark.parametrize("g,x,h,p", bsgs)
def test_baby_step_giant_step(g, x, h, p):
    x = discrete_log.baby_step_giant_step(g, h, p, p - 1)
    assert x != 0 and x is not None
    assert g**x % p == h

def test_pohlig_hellman():
    n = 101
    f = {2: 2, 5: 2}
    g, h = 15, 100
    assert discrete_log.pohlig_hellman(g, h, n, f) == 50


def test_pohlig_hellman_prime_power():
    g, h, p, e, o = 27, 40, 2, 3, 41
    assert discrete_log.pohlig_hellman_prime_power_order(g, h, p, e, o) == 4
