from mobius import GeneralizedContinuedFraction, MobiusTransform
from functools import reduce
from math import factorial
import numpy as np
import operator
import mpmath
import sympy
from convergence_rate import calculate_convergence


def calculate_lhs_family1_non_matrix(k):
    def calculate_d():
        d_ = [0, 1]
        for n in range(1, k + 1):
            d_.append((8 * n ** 2 + 2 * n + 1) * d_[n] - 2 * n * (2 * n - 1) ** 3 * d_[n - 1])
        return -d_[k - 1]

    def calculate_d_closed():
        with mpmath.workdps(10):
            _a = (mpmath.gamma(k - 0.5) ** 2) / (4 * mpmath.pi)
            _b = mpmath.power(2, 2 * k + 1) * mpmath.catalan
            _c = mpmath.pi * mpmath.gamma(2 * k - 1) * (
                    mpmath.hyper([1, k - 0.5, k], [k + 0.5, k + 0.5], 1) / (mpmath.gamma(k + 0.5) ** 2))
        return -_a * (_b - _c)

    b = factorial(2 * k - 2)
    c = 2 * reduce(operator.mul, [i for i in range(-1, 2 * k - 1, 2)]) ** 2
    d = calculate_d()
    # d = mpmath.nint(calculate_d_closed())
    lhs = b / (c * mpmath.catalan + d)
    lhs_sym = b / (c * sympy.Catalan + d)
    return lhs, lhs_sym


def calculate_lhs_family1(k):
    k = k - 1  # for calculation
    func = MobiusTransform(np.array([[2, 0], [0, 1]]))
    for n in range(1, k + 1):
        func = MobiusTransform(np.array([[(2 * n - 1) ** 2, -1], [0, 2 * n * (2 * n - 1)]], dtype=int)) * func
    func = func.reciprocal()
    return func(mpmath.catalan), sympy.sympify(func.sym_expression(sympy.Catalan))


def calculate_lhs_family2(k):
    k = k - 1
    func = MobiusTransform(np.array([[2, -1], [0, 2]]))
    for n in range(1, k+1):
        func = MobiusTransform(np.array([[(2 * n - 1) * (2 * n - 3), -1], [0, 2 * n * (2 * n - 3)]], dtype=int)) * func
    func = func.reciprocal()
    return func(mpmath.catalan), sympy.sympify(func.sym_expression(sympy.Catalan))


def check_catalan(an_func, bn_func, calculate_lhs):
    for N in range(1, N_to_check + 1):
        an = [an_func(n, N) for n in range(1000 + N * 100)]
        bn = [bn_func(n, N) for n in range(1000 + N * 100)]
        with mpmath.workdps(400):
            lhs, lhs_sym = calculate_lhs(N)
            gcf = GeneralizedContinuedFraction(an, bn)
            rhs_val = gcf.evaluate()

            lhs_str = mpmath.nstr(lhs, 100)
            rhs_str = mpmath.nstr(rhs_val, 100)
            eq = sympy.Eq(lhs_sym, gcf.sym_expression(4))
            sympy.pprint(eq)
            if lhs_str == rhs_str:
                print('convergence rate: {}'.format(mpmath.nstr(calculate_convergence(gcf, lhs), 5)))
                print('equal!\n')
            else:
                print('error\n')


N_to_check = 8


def family1_check():
    # family 1 - (3,3,1), (-2,1,1,1,N)
    check_catalan(lambda n, N: n * (3 * n + 3) + 1,
                  lambda n, N: -2 * ((n + 1) ** 3) * (n + N),
                  calculate_lhs_family1)


def family1_check_non_matrix():
    # family 1 - (3,3,1), (-2,1,1,1,N) - non matrix lhs form
    check_catalan(lambda n, N: n * (3 * n + 3) + 1,
                  lambda n, N: -2 * ((n + 1) ** 3) * (n + N),
                  calculate_lhs_family1_non_matrix)


def family2_check():
    # family 2 - (3,7,3), (-2,1,1,3,N)
    check_catalan(lambda n, N: n * (3 * n + 7) + 3,
                  lambda n, N: -2 * ((n + 1) ** 2 * (n + 3)) * (n + N),
                  calculate_lhs_family2)


# family1_check()
# family1_check_non_matrix()
family2_check()
