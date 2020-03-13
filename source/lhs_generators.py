from time import time
import pickle
import itertools
from sympy import Abs, symbols
from enumerate_over_signed_rcf import SignedRcfEnumeration

"""
The purpose of this file is to allow for the construction enumerations following a specific pattern. 
Similar to the series generators of the MITM flavor. Here we chose patterns for the LHS.
The enumeration of the LHS is a substantial part of the search. 
"""

def create_standard_lhs(poly_deg, coefficients_limit, out_path=None, do_print=True):
    """
    Prepares generic LHS enumerations (variable 'x'). When used later in searching conjectures- 'x' will be substituted
    to any constant.
    :param poly_deg: Degree of polynomials in rational function.
    :param coefficients_limit: Limit for coefficients (symmetrical).
    :param out_dir: Directory for saving the binary LHS enumeration.
    :return: List of standard generic LHS enumerations that can later be used for enumerating over any constant.
    """
    if do_print:
        print("starting enumeration")
    strt = time()
    x = symbols('x')
    enum = SignedRcfEnumeration(sym_constant=x, cycle_len_range=None, coefficients_limit=coefficients_limit,
                                poly_deg=poly_deg, do_print=do_print)
    generic_variations = enum.create_rational_variations_enum()
    if out_path is not None:
        with open(out_path, 'wb') as file:
            pickle.dump(generic_variations, file)
    if do_print:
        print("Finished. Took {} sec".format(time() - strt))
    return generic_variations


def create_biased_monoms(max_deg, const_coeff_lim, bias_lim):  # monom = x**d.
    """
    Prepares a generic LHS enumeration of the form: ax**k+b. Where a,b,k  are integers, and x is the variable.
    :param max_deg: maximum degree (and minimum since symmetrical) of the monom.
    :param const_coeff_lim: symmetrical limit for the coefficient of the power.
    :param bias_lim: symmetrical limit for the coefficient of the bias.
    :return: List of custom generic LHS enumerations that can later be used for enumerating over any constant.
    """
    x = symbols('x')
    c_list = [j for j in range(1, const_coeff_lim + 1) if j != 0]
    d_list = [j for j in range(-max_deg, max_deg + 1) if j != 0]
    b_list = [j for j in range(-bias_lim, bias_lim + 1)]
    domain = itertools.product(c_list, d_list, b_list)
    biased_monoms = []
    for it in domain:
        biased_monoms.append(Abs(it[0] * x ** it[1] + it[2]))
    return biased_monoms


def e_to_the_k_hypo(max_deg, bias_lim):
    x = symbols('x')
    d_list = [j for j in range(-max_deg, max_deg + 1) if j != 0]
    b_list = [j for j in range(-bias_lim, bias_lim + 1)]
    domain = itertools.product(d_list, b_list)
    variations = []
    for it in domain:
        variations.append(((Abs(it[0]) - 1)*(x**it[0])) + it[1])
    return variations


if __name__ == '__main__':
    print("Keep it simple, use the API")
