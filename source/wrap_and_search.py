from massey import create_series_from_shift_reg, slow_massey

import itertools
import math
import sympy
import mpmath
from sympy import lambdify
from mpmath import mpf as dec
from ortools.linear_solver.pywraplp import Solver
from mobius import MobiusTransform, SimpleContinuedFraction, GeneralizedContinuedFraction
import numpy as np
from sympy import pi
from sympy import E as e
from sympy import besseli
from sympy import lambdify
from sympy import zeta
import sympy
from operator import itemgetter


# flavors of find transform:
def find_transform_2const_lin(xx1, xx2, y, limit, threshold=1e-7):
    """
    find a integer solution to a1x1 + a2x2 + b - c1x1y - c2x2y - dy = 0
    this will give us the mobius transform: T(x) = y
    :param x: numeric constant to check
    :param y: numeric manipulation of constant
    :param limit: range to look at
    :return MobiusTransform in case of success or None.
    """
    x1 = xx1
    x2 = xx2
    x3 = dec(1.0)
    x4 = -xx1*y
    x5 = -xx2*y
    x6 = -y
    solver = Solver('mobius', Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    a1 = solver.IntVar(-limit, limit, 'a1')
    a2 = solver.IntVar(-limit, limit, 'a2')
    b = solver.IntVar(-limit, limit, 'b')
    c1 = solver.IntVar(-limit, limit, 'c1')
    c2 = solver.IntVar(-limit, limit, 'c2')
    d = solver.IntVar(-limit, limit, 'd')
    f = solver.NumVar(0, 1, 'f')
    solver.Add(f == (a1*x1 + a2*x2 + b*x3 + c1*x4 + c2*x5 + d*x6))
    solver.Add(a1*x1 + a2*x2 + b >= 1)   # don't except trivial solutions and remove some redundancy
    solver.Minimize(f)
    status = solver.Solve()
    if status == Solver.OPTIMAL:
        if abs(solver.Objective().Value()) <= threshold:
            return int(a1.solution_value()), int(a2.solution_value()), \
                   int(b.solution_value()), int(c1.solution_value()), \
                   int(c2.solution_value()), int(d.solution_value())
    else:
        return None


def is_good_gcf(a_, b_, integer_limit, constants):
    """
    Finds conjectures between a GCF and mobius transforms of math constants.
    :param a_: sequence
    :param b_: sequence
    :param integer_limit: maximum allowed value for mobius integers.
    :param constants: list of math fundemental constants.
    :return: polynomial coefficients of P field.
    """
    with mpmath.workdps(1000):
        y = GeneralizedContinuedFraction(a_, b_).evaluate()
        good_results = []
        for x in constants:
            xx = lambdify((), x, modules="mpmath")()
            yy = lambdify((), y, modules="mpmath")()
            result = find_transform_2const_lin(xx, yy, integer_limit)
            if result is not None:
                good_results.append([a_, b_, x, result])
        return good_results



"""
Follows a class for enumeration of signed regular CFs. 
There appears to be some redundancy with signed RCFs fitting the same lhs. will look into it more.
Further documentation to follow. 
unittest also added to tests. 
"""
class SignedRcfEnumeration(object):
    def __init__(self, sym_constant, coefficients_limit, depth, prime=199):
        self.beauty_standard = 30
        self.enum_dps = 50
        self.verify_dps = 200
        self.coeff_lim = coefficients_limit
        self.const_sym = sym_constant
        self.const_val = lambdify((), sym_constant, modules="mpmath")
        self.depth = depth
        self.prime = prime

    def create_sign_seq_enumeration(self, cycle_len):
        return [list(sign_seq) for sign_seq in list(itertools.product([-1,1], repeat=cycle_len))]

    def create_rational_variations_enum(self, poly_deg):
        coeffs = [i for i in range(-self.coeff_lim, self.coeff_lim+1)]
        numerators = [list(numer) for numer in list(itertools.product(coeffs, repeat=poly_deg+1))]
        denominators = [list(denom) for denom in list(itertools.product(coeffs, repeat=poly_deg+1))]
        return [list(var) for var in list(itertools.product(numerators,denominators))]

    def find_signed_rcf_conj(self, poly_deg, max_cycle_len):
        results = []
        rational_variations = self.create_rational_variations_enum(poly_deg)
        for seq_len in range(2, max_cycle_len+1):
            sign_cycles = self.create_sign_seq_enumeration(seq_len)
            for rat_var in rational_variations:
                numer, denom = [], []
                for i in range(len(rat_var[0])):
                    numer.append(rat_var[0][i])
                    denom.append(rat_var[1][i])
                pad_len = 6-len(numer)
                pad = [0]*pad_len
                numer = numer + pad
                denom = denom + pad
                if numer == denom or denom == [0, 0, 0, 0, 0, 0]:
                    continue

                var_sym = (numer[0] + numer[1]*self.const_sym + numer[2]*self.const_sym**2 +
                           numer[2]*self.const_sym**3 + numer[4]*self.const_sym**4 +
                           numer[5]*self.const_sym**5) / \
                          (denom[0] + denom[1]*self.const_sym + denom[2]*self.const_sym**2 + denom[3]*self.const_sym**3
                           + denom[4]*self.const_sym**4 + denom[5]*self.const_sym**5)

                var_gen = lambdify((), var_sym, modules="mpmath")
                for sign_cyc in sign_cycles:
                    temp = 2
                    new_cyc = True
                    while temp < len(sign_cyc) and new_cyc:
                        if len(sign_cyc) % temp == 0:
                            short_cyc = sign_cyc[temp:]
                            if sign_cyc == short_cyc * (len(sign_cyc) // temp):
                                new_cyc = False
                            temp *= 2
                        else:
                            temp *= 2
                    if not new_cyc:
                        continue

                    b_ = (sign_cyc * (self.depth//seq_len))[:self.depth]
                    signed_rcf = GeneralizedContinuedFraction.from_irrational_constant(const_gen=var_gen, b_=b_)
                    a_ = signed_rcf.a_
                    if 0 in a_:
                        continue
                    if len(a_) < self.depth:
                        continue
                    leo = slow_massey(a_, self.prime)
                    if len(leo) < self.beauty_standard:
                        results.append([numer, denom, sign_cyc, a_])
        return results


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
    cart_prod = list(itertools.product(rng, repeat=reg_size))
    polynoms = [list(prod) for prod in cart_prod]
    registers = [list(prod) for prod in cart_prod]
    for poly in polynoms:
        poly.insert(0, 1)
    # now we have a list of polynoms, and a list of registers(start conditions)
    # for each iteration we want a_pol, b_pol, a_reg, b_reg.
    domain = list(itertools.product(polynoms, polynoms, registers, registers))  # list of tuples
    # every element in domain is 2 polynoms, and 2 start registers. Starting search:
    possible_conjectures = []
    for instance in domain:
        [poly_a, poly_b, reg_a, reg_b] = list(instance)
        results = is_good_gcf(create_series_from_shift_reg(poly_a, reg_a, seq_depth),
                             create_series_from_shift_reg(poly_b, reg_b, seq_depth), int_lim, constants)
        # check if results from this point are massey-pretty:
        if results is not None:
            for res in results:
                a_sr = slow_massey(res[0], prime)
                b_sr = slow_massey(res[1], prime)
                if len(a_sr) <= thresh and len(b_sr) <= thresh:
                    possible_conjectures.append(res)

#
# if __name__ == '__main__':
#     enumerator = SignedRcfEnumeration(e, 1, 100, 199)
#     results = enumerator.find_signed_rcf_conjs(2, 2)
#     print("found {} results for e:")
#     for res in results: