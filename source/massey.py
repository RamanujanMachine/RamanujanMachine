# Dan's advice on choosing the prime:
# Regarding choice of primes - I dunno :) something large.
# Things to worry about are the prime being smaller than K+1,
# [[K is the maximum order of a polynomial on one of the channels]]
# or dividing the leading coefficients of some of the P_m     [[a_n = a_(N*m + j) = P_j(m)]]
# (or perhaps a combination of thereof), and therefore decreasing the degree.


"""
next step:
    make it fast. how?
        - convert it to Class, keep table of inverts
        - avoid unnecessary allocations. keep max poly degree?
        - vectorize?
        - faster modulo?
TODO - unittests
"""

from numpy import array, int64, concatenate


def _inv_mod(a, p):  # get inverted modulo of prime-field p.
    b = p
    caa, cba = 1, 0
    while b != 0:
        caa, cba = cba, caa - (a // b) * cba
        a, b = b, a % b
    return caa


def _poly_add(a_x, b_x):  # add 2 polynomials
    if len(a_x) > len(b_x):
        b_x = concatenate((b_x, [0] * (len(a_x) - len(b_x))))
    elif len(a_x) < len(b_x):
        a_x = concatenate((a_x, [0] * (len(b_x) - len(a_x))))
    return a_x + b_x


def _update_polynomial(c_x, b_x, d, b, m, p):  # update step of massey algorithm
    b_x = concatenate(([0] * m, b_x))
    q = (d * _inv_mod(b, p)) % p
    b_x = (-b_x * q) % p
    c_x = _poly_add(c_x, b_x) % p
    return c_x


def slow_massey(line, p):
    """
    Apply "Berlekamp-Massey" Algorithm on series.
    :param line: input series
    :param p: prime number field
    :return: polynomial coefficients of P field.
    """
    assert p < 2 ** 32
    s_ = array(line, dtype=int64)  # input series
    c_ = array([1], dtype=int64)  # current polynomial
    b_ = array([1], dtype=int64)  # previous error polynomial
    poly_deg = 0  # current polynomial degree
    m = 1  # number of iterations since last error
    b = 1  # copy of the last discrepancy d

    for n in range(0, len(line)):
        d = s_[n]
        for i in range(1, poly_deg + 1):
            d += (c_[i] * s_[n - i]) % p
            d %= p
        if d == 0:
            m += 1
        elif (2 * poly_deg) <= n:
            temp_copy = c_.copy()
            c_ = _update_polynomial(c_, b_, d, b, m, p)
            poly_deg = n + 1 - poly_deg
            b_ = temp_copy.copy()
            b = d
            m = 1
        else:
            c_ = _update_polynomial(c_, b_, d, b, m, p)
            m = m + 1
    try:
        first_non_zero = [x for x in c_ if x != 0][0]
        c_ = (c_ * _inv_mod(first_non_zero, p) + (p // 2)) % p - (p // 2)
    except RuntimeError:  # Exception will only happen if c_ is all zeros, which is bad anyway
        print('Something went wrong with slow_massey algorithm')
    return c_


def create_series_from_polynomial(poly_a, n):
    """
    given a polynomial P(x) as generator return An when A[i] = P(i)
    :param poly_a:
    :type poly_a: list
    :param n: number of values to generate
    :type n: int
    :return: series with length of n
    """
    a_ = []
    for i in range(n):
        a_i = 0
        pol = 1
        for j in range(len(poly_a)):
            a_i += pol * poly_a[j]
            pol *= i
        a_.append(a_i)
    return a_


def create_series_from_shift_reg(poly_a, initials, n):
    """
    given a polynomial P as shift register, generate a series.
    series follows the following rule: An = sum 1 to deg(P): (P(i) * A(n-i))
    :param poly_a: P(x) : same format as the massey algorithm output (expects poly_a[0] = 1)
    :param initials: starting conditions of the series. must be of length = deg(P)
    :param n: number of values to generate
    :return: series with length of n
    """
    assert ((len(poly_a) - len(initials)) == 1), 'There should be deg(poly) initial conditions'
    assert (len(poly_a) > 0) and (poly_a[0] == 1), 'Illegal polynomial - first coefficient must be 1'
    a_ = []
    for i in range(len(initials)):
        a_.append(initials[i])
    for i in range(len(a_), n):
        a_i = 0
        for j in range(1, len(poly_a)):
            a_i -= poly_a[j] * a_[i - j]
        a_.append(a_i)
    return a_
