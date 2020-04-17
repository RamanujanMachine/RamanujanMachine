import numpy as np
from typing import List


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


# an_coefs = [n * (n * (6 * n + 9) + 5) + 1 for n in range(10)]
# print(find_polynomial_series_coefficients(3, an_coefs, 0))
