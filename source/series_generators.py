import sympy
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
    this is the reversed action to the massey algorithm.
    given a polynomial P as shift register, generate a series.
    series follows the following rule: A[n] = sum i=1 to deg(P): -(P[i] * A[n-i])
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


def create_series_from_compact_poly(poly_a, n):
    """
    create a series of type n(n(...(a[0]*n + a[1]) + a[2]) + ...) + a[k]
    :param poly_a: a[k] coefficients
    :param n: length of series
    :return: a list of numbers in series
    """
    ret = []
    for i in range(n):
        tmp = 0
        for c in poly_a:
            tmp *= i
            tmp += c
        ret.append(tmp)
    return ret


def create_series_from_compact_poly_with_shift1(poly_a, n):
    """
    create a series of type m(m(...(a[1]*m + a[0]) + a[2]) + ...) + a[k],
    where m=n+1
    :param poly_a: a[k] coefficients
    :param n: length of series
    :return: a list of numbers in series
    """
    ret = []
    for m in range(1, n+1):
        tmp = 0
        for c in poly_a:
            tmp *= m
            tmp += c
        ret.append(tmp)
    return ret


def create_zeta_bn_series(deg, poly_a, n):
    """
    create a series of type: a[0]*(n+1)^d - a[1]*(n+1)^(d-1)
    apparently, this is useful for building zeta function CFs
    :param deg: d
    :param poly_a: coefficients a[0], a[1]. if more exist, they are ignored.
    :param n: depth of series
    :return: a list of numbers in series
    """
    ret = []
    for i in range(n):
        res = poly_a[0]*((i+1)**deg) - poly_a[1]*((i+1)**(deg-1))
        ret.append(res)
    return ret

def zeta3_an_generator(x_, n):
    """
    Generator3[x3_, x0_] := {x0, 2 *x0 + x3, 3*x3, 2*x3}
    :param x_:
    :param n:
    """
    ret = []
    for i in range(n):
        res = x_[0] + (2*x_[0] + x_[1])*i + 3*x_[1]*(i**2) + 2*x_[1]*(i**3)
        ret.append(res)
    return ret


def zeta5_an_generator(x_, n):
    """
    Generator5[x5_, x3_, x0_] := {x0, 2 *x0 + x3 - 2 *x5, 3*x3 - 5*x5, 2*x3, 5*x5, 2*x5}
    :param x_:
    :param n:
    """
    ret = []
    for i in range(n):
        res = x_[0] + (2*x_[0] + x_[1] - 2*x_[2])*i + (3*x_[1] - 5*x_[2])*(i**2) +\
              2*x_[1]*(i**3) + 5*x_[2]*(i**4) + 2*x_[2]*(i**5)
        ret.append(res)
    return ret


def catalan_bn_generator(x_, n):
    ret = []
    for i in range(n):
        res = x_[0]*(2*i+1)**4 + x_[1]*(2*i+1)**3
        ret.append(res)
    return ret