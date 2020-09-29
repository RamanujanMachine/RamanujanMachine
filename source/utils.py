import numpy as np
from typing import List
import time
import mpmath
import matplotlib.pyplot as plt

# Measures the amount of time the function takes to run in milliseconds in order to check improvements
def measure_performance(func):
    """
    measure the duration of a function execution in milliseconds in order to improve its performance

    Use as decorator by adding @measure_performance before functions you wish to test
    """
    def wrapper(*args, **kwargs):
        start = time.time()
        func_value = func(*args, **kwargs)
        end = time.time()
        print(f"{func.__name__} elapsed (with compilation) = {1000 * (end - start)} ms")
        return func_value

    return wrapper


def find_polynomial_series_coefficients(poly_deg, lead_terms: List[int], starting_n=0):
    """
    find polynomial coefficients of polynomial integer series, by solving linear equations.
    :param poly_deg: degree of polynomial
    :param lead_terms: list of the first few terms of the series (must be at least poly_deg+1)
    :param starting_n: what index does the series start from (default is 0).
    :return: list of polynomial coefficients (starting with the lead term)
    """
    assert (poly_deg+1) <= len(lead_terms)
    mat = np.array([[n**i for i in range(poly_deg+1)] for n in range(starting_n, starting_n+poly_deg+1)])
    x = np.array(lead_terms[:poly_deg+1])
    coefs = np.linalg.inv(mat) @ x
    coefs = np.flip(coefs)
    int_coefs = np.round(coefs)
    if any(np.abs(int_coefs-coefs) > 10e-7):
        print('warning! non integer coefficients - {}'.format(coefs))
        return list(coefs)
    else:
        return list(int_coefs.astype(int))


def get_poly_deg_and_leading_coef(poly_coef):
    deg = len(poly_coef) - 1
    for i in poly_coef:
        if i != 0:
            return deg, i
        else:
            deg -= 1
    # if we get here the degree is 0
    return deg, poly_coef[-1]


def iter_series_items_from_compact_poly(poly_a, max_runs, starting_n=0):
    """
    create a series of type n(n(...(a[0]*n + a[1]) + a[2]) + ...) + a[k]
    :param poly_a: a[k] coefficients
    :param n: length of series
    :return: a list of numbers in series
    """
    for i in range(starting_n, max_runs + starting_n):
        tmp = 0
        for c in poly_a:
            tmp *= i
            tmp += c
        yield tmp


def plot_gcf_convergens(an_poly_coef, bn_poly_coef, max_iters, divide_interval=101, label=None):
    computed_values = []
    label = f'an {an_poly_coef} bn {bn_poly_coef}' if not label else label

    # a_n iterator is always one ahead, while a[0] goes to p. so this is 
    # an ugly hack

    poly_a_deg, _ = get_poly_deg_and_leading_coef(an_poly_coef)
    poly_b_deg, _ = get_poly_deg_and_leading_coef(bn_poly_coef)

    print(f'deg a - {poly_a_deg} deg b - {poly_b_deg}')
    if poly_a_deg * 2 == poly_b_deg:
        print('\texpo')
        print(f'\t4b_d = {bn_poly_coef[0] * 4} | -a^2 = {-1 * (an_poly_coef[0]**2)}')
        if bn_poly_coef[0] * 4 > -1 * (an_poly_coef[0]**2):
            print('\t\tcond passed')
        elif bn_poly_coef[0] * 4 == -1 * (an_poly_coef[0]**2):
            print('\t\tequality')
        else:
            print('\t\tcond failed')
    elif poly_a_deg * 2 > poly_b_deg:
        print('\tsuper expo')
    else:
        print('\tsub expo')
    an_items_iterator = iter_series_items_from_compact_poly(an_poly_coef, max_iters)
    bn_items_iterator = iter_series_items_from_compact_poly(bn_poly_coef, max_iters, 1)
    
    prev_q = 0
    q = 1
    prev_p = 1
    p = an_items_iterator.__next__() # will place a[0] to p

    for i, (a_i_1, b_i) in enumerate(zip(an_items_iterator, bn_items_iterator)):
        # a_i_1 is the (i+1)'th item of an, and b_i the the i'th item of bn
        tmp_a = q
        tmp_b = p

        q = a_i_1 * q + b_i * prev_q
        p = a_i_1 * p + b_i * prev_p
        
        prev_q = tmp_a
        prev_p = tmp_b

        # This is the hard part to compute.
        # TODO - consider computing once in every X iters; Might be significant 
        if i % divide_interval != 0:
            continue
        if q != 0:  # safety check
            computed_values.append((mpmath.mpf(p) / mpmath.mpf(q), i))
        else:
            computed_values.append((mpmath.mpf(0), i))

    global current_fig
    current_fig += 1 
    plt.figure(current_fig)
    plt.plot([i[1] for i in computed_values], [i[0] for i in computed_values], '+')
    plt.ion()
    plt.title(label=label)
    plt.show()

    return computed_values
