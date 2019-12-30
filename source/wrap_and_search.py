import mobius
from massey import create_series_from_shift_reg, slow_massey
import itertools
import math
import sympy
import mpmath
from sympy import lambdify
from mobius import is_good_gcf

def find_gcf_conj(gen_pol_lim, int_lim, constants, reg_size, seq_depth):
    """
      Looks for conjecture whose sequences are massey-pretty.
      :param gen_pol_lim: maximum coefficient allowed for the polynomials generating the sequences.
      :param int_lim: maximum integer allowed for mobius transforms in conjectures.
      :param reg_size: length of shift register generating sequences. polynomial length always 1 above.
      :param constants: list of math fundemental constants.
      :param seq_depth: depth of GCF to be developed before checking results
      :return: list of conjectures in the form: [a_, b_, constant, mobius]
      """
    prime = 199
    thresh = 30

    rng = range(gen_pol_lim)
    cart_prod = list(itertools.product(rng, repeat=reg_size))  # list of tuples
    polynoms = [list(prod) for prod in cart_prod]  # lists of lists
    registers = [list(prod) for prod in cart_prod]  # lists of lists
    for poly in polynoms:
        poly.insert(0, 1)
    # now we have a list of polynoms, and a list of registers(start conditions)
    # for each iteration we want a_pol, b_pol, a_reg, b_reg.
    domain = list(itertools.product(polynoms, polynoms, registers, registers))  # list of tuples
    # every element in domain is 2 polynoms, and 2 start registers. Starting search:
    possible_conjectures = []
    for instance in domain:
        [poly_a, poly_b, reg_a, reg_b] = list(instance)  # list of lists.
        results = is_good_gcf(create_series_from_shift_reg(poly_a, reg_a, seq_depth),
                             create_series_from_shift_reg(poly_b, reg_b, seq_depth), int_lim, constants)
        # check if results from this point are massey-pretty:
        if results is not None:
            for res in results:
                #prime = get_prime_greater(max(max(res[0]), max(res[1])))
                a_sr = slow_massey(res[0], prime)
                b_sr = slow_massey(res[1], prime)
                if len(a_sr) <= thresh and len(b_sr) <= thresh:
                    possible_conjectures.append(res)
