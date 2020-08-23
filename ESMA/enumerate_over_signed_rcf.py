import os
import pickle
from time import time
import itertools
import mpmath
import sympy
from sympy import lambdify, Rational
from massey import slow_massey
from EfficientGCF import EfficientGCF
from ramanujan.utils.mobius import GeneralizedContinuedFraction
from ramanujan.utils.convergence_rate import calculate_convergence

"""
Some important terminology:
LFSR - Linear Feedback Shift Register. With initial conditions, it can represent a variety of sequences over any field.
       Form of the  output of Berlekamp–Massey algorithm. A finite LFSR indicates sequence has an underlying structure.
LHS - Some parametric function f(constant, coefficients). By default, rational functions with integer coefficients.
      However, any parametric function can be used by writing a fitting generator (See lhs_generators.py)
"""

def clear_end_zeros(items):
    """
    removes zeros at end of list.
    This is sometime required when Massey output has zeros at the end.
    :param items: iterable.
    :return: list without the zeros in the end.
    """
    while items[-1] == 0:
        items.pop()


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


class SignedRcfEnumeration(object):

    def __init__(self, sym_constant, cycle_len_range, depth=100, coefficients_limit=None, poly_deg=None, min_deg=None,
                 prime=199, custom_enum=None, do_print=True):
        """
        Initialize search engine.
        Basically, this is a 3 step procedure:
        1) Enumerates LHS symbolic expressions of functions of the constant, and non repeating sign periods.
        2) Iterates through domain. With low precision extracts a series. Checks if massey-pretty. Saves hits.
        3) Refine results - takes results from (2) and validate them to 100 decimal digits.
        Note that the structure of the enumeraion (and in fact of the problem), makes it possible to divide
        any domain to separate domains with respect to the signed periods, polynom degrees, but not the coefficients
        limit.
        Custom parametrized functions can be used. In that case the LHS enumeration must be generated externally
        (see lhs_generators.py).
        :param sym_constant: sympy constant
        :param cycle_len_range: range of lengths for the sign sequence's period.
        :param depth: Number of elements of a series to extract. Relates to length of typical LFSRs  of the consant.
        :param coefficients_limit: Range of coefficients for the rational function on the LHS.
        :param poly_deg: Maximum degree of numerator and denominator polynomials in the rational LHS.
        :param min_deg: Used to exclude lower degree numerator and denominator polynomials in rational LHS from search.
        :param prime: Prime number in use by Massey algorithm.
        :param custom_enum: A ready-made enumeration that only requires substituting a variable 'x' with the constant.
        :param do_print: Print outputs (Used as False primarily for unit tests).
        """
        self.enum_dps = 500
        self.verify_dps = 1000
        self.coeff_lim = coefficients_limit
        if cycle_len_range is not None:
            self.min_cycle_len = cycle_len_range[0]
            self.max_cycle_len = cycle_len_range[1]
        self.poly_deg = poly_deg
        self.min_deg = min_deg
        self.const_sym = sym_constant
        try:
            self.const_val = lambdify((), sym_constant, modules="mpmath")
        except AttributeError:
            self.const_val = self.const_sym.mpf_val

        self.depth = depth
        self.beauty_standard = self.depth//4
        """
        'beauty standard' is a Threshold for the length of a sequence LFSR (Berlekamp–Massey output). Only extracted 
        sequences with shorter LFSR will be considered as a possible valuable result.
        """
        self.verify_depth = 1000
        self.prime = prime
        self.custom_enum = custom_enum
        self.do_print = do_print

    def create_sign_seq_enumeration(self):
        """
        Creates a list of all possible sign sequences.
        Uses either a pre defined dictionary or a function to skip redundancies.
        """
        sign_seqs = []
        for cyc_len in range(self.min_cycle_len, self.max_cycle_len+1):
            sign_seqs = sign_seqs + list(itertools.product([-1, 1], repeat=cyc_len))
        return sign_seqs

    def create_rational_symbol(self, numerator, denominator):
        """
        creates a symbolic expression of a rational expression: P(c)/Q(c). Where P,Q polynomials, and c is the constant.
        :param numerator: numerator polynomial coefficients where numerator[0] is the free coefficient.
        :param denominator: denominator polynomial coefficients where denominator[0] is the free coefficient.
        :return: sympy symbolic expression.
        """
        numer_deg = len(numerator) - 1
        denom_deg = len(denominator) - 1
        numer_sym = 0
        denom_sym = 0
        for i in range(numer_deg + 1):
            numer_sym += numerator[i]*(self.const_sym**i)
        for i in range(denom_deg + 1):
            denom_sym += denominator[i]*(self.const_sym**i)
        return numer_sym/denom_sym

    def create_rational_variations_enum(self):
        """
        Creates a list of all possible rational expressions for the LHS.
        Expressions saved as sympy-simplified, positive expressions to reduce redundancy.
        Additional checks are performed to exclude degenerated cases.
        """
        if self.do_print:
            print("Starting enumeration over LHS")
        unsimplified = (2*self.coeff_lim + 1) ** (2*(self.poly_deg + 1))
        if self.do_print:
            print("Number of variations before simplification is: {}".format(unsimplified))
        start = time()
        coeffs = [i for i in range(-self.coeff_lim, self.coeff_lim+1)]
        if self.min_deg is not None:
            numerators = [list(numer) for numer in list(itertools.product(coeffs, repeat=self.poly_deg + 1))
                          if len(numer) >= self.min_deg+1]
            denominators = [list(denom) for denom in list(itertools.product(coeffs, repeat=self.poly_deg + 1))
                            if len(denom) >= self.min_deg+1]
        else:
            numerators = [list(numer) for numer in list(itertools.product(coeffs, repeat=self.poly_deg+1))]
            denominators = [list(denom) for denom in list(itertools.product(coeffs, repeat=self.poly_deg+1))]
        variations = itertools.product(numerators, denominators)
        expressions = set()
        cnt = 0
        mytimer= time()
        for var in variations:
            if ((cnt % 1000) == 0) and (cnt != 0) and self.do_print:
                print("{} variations took {} minutes".format(cnt, round((time()-mytimer)/60, 2)))
            cnt += 1
            numer = var[0]
            denom = var[1]
            if (denom == [0 for i in denom]) or (numer == [0 for i in numer]):
                continue
            var_sym = self.create_rational_symbol(numer, denom)
            var_sym = sympy.simplify(var_sym)
            if isinstance(var_sym, Rational):  # Simplified expression is a rational number (independant of constant).
                continue
            if abs(var_sym) not in expressions:
                expressions.add(abs(var_sym))
        if self.do_print:
            ("Finished enumerations. Took {}  seconds".format(round(time()-start, 2)))
        return expressions

    def find_signed_rcf_conj(self):
        """
        Builds the final domain.
        Iterates throgh the domain:
        extraction->massey->check->save.
        Additional checks are performed to exclude degenerated cases.
        If a generic enumeration is given will use it instead of enumerating.
        """
        inter_results = []
        redundant_cycles = set()
        # Enumerate:
        if self.custom_enum is None:
            lhs = self.create_rational_variations_enum()
        else:
            if self.do_print:
                print("Substituting " + str(self.const_sym) + ' into generic LHS:')
            strt = time()
            lhs = [var.subs({sympy.symbols('x'): self.const_sym}) for var in self.custom_enum]
            if self.do_print:
                print("Took {} sec".format(time() - strt))
        sign_seqs = []
        for cyc_len in range(self.min_cycle_len, self.max_cycle_len + 1):
            sign_seqs = sign_seqs + list(itertools.product([-1, 1], repeat=cyc_len))
        domain_size = len(lhs) * len(sign_seqs)
        if self.do_print:
            print("De-Facto Domain Size is: {}\n Starting preliminary search...".format(domain_size))
        checkpoint = max(domain_size // 20, 5)
        count = 0
        start = time()
        # Iterate
        bad_variation = []
        for instance in itertools.product(lhs, sign_seqs):
            count += 1
            var, sign_period = instance[0], list(instance[1])
            if var == bad_variation:
                continue
            bad_variation = []
            if ''.join([str(c) for c in sign_period]) in redundant_cycles:
                continue
            # if this cycle was not redundant it renders some future cycles redundant:
            for i in range(2, (self.max_cycle_len // len(sign_period)) + 1):
                redun = sign_period * i
                redundant_cycles.add(''.join([str(c) for c in redun]))
            var_gen = lambdify((), var, modules="mpmath")
            seq_len = len(sign_period)
            if (count % checkpoint == 0) and (self.do_print):
                print("\n{}% of domain searched.".format(round(100 * count / domain_size, 2)))
                print("{} possible results found".format(len(inter_results)))
                print("{} minutes passed.\n".format(round((time() - start) / 60, 2)))
            b_ = (sign_period * ((self.depth // seq_len) + 1))  # Concatenate periods to form sequence.
            b_ = b_[:self.depth]  # Cut to proper size.
            with mpmath.workdps(self.enum_dps):
                try:
                    signed_rcf = GeneralizedContinuedFraction.from_irrational_constant(const_gen=var_gen, b_=b_)
                except ZeroDivisionError:
                    if self.do_print:
                        print('lhs:')
                    sympy.pprint(var)
                    bad_variation = var
                    continue
            a_ = signed_rcf.a_
            if 0 in a_:
                continue
            if len(a_) < self.depth:
                continue
            a_lfsr = list(slow_massey(a_, self.prime))
            clear_end_zeros(a_lfsr)
            if len(a_lfsr) < self.beauty_standard:
                inter_results.append([var, sign_period, a_[:(len(a_lfsr)-1)], a_lfsr])
        return inter_results

    def verify_results(self, results):
        """
        Validate intermediate results to 100 digit precision
        If a numeric value appears multiple times, the first is kept as valid. The rest saved as recurring for later.
        """
        verified_results = []
        recurring_value_results = {}
        res_set = set()
        for res in results:
            var_gen = lambdify((), res[0], modules="mpmath")
            a_ = create_series_from_shift_reg(res[3], res[2], self.verify_depth)
            b_ = (res[1] * ((self.verify_depth // len(res[1])) + 1))
            b_ = b_[:self.verify_depth]
            gcf = EfficientGCF(a_, b_)
            with mpmath.workdps(self.verify_dps):
                lhs_str = mpmath.nstr(var_gen(), 100)
                rhs_val = gcf.evaluate()
                rhs_str = mpmath.nstr(rhs_val, 100)
                if rhs_str != lhs_str:
                    continue
                key = lhs_str
            if key not in res_set:
                res_set.add(key)
                verified_results.append(res)
                recurring_value_results[key] = []
            else:
                recurring_value_results[key].append(res)
                continue
        return verified_results, recurring_value_results

    def print_results(self, results, latex=True):
        """
        Print results in either unicode or LaTex.
        :param results: verified results.
        :param latex: LaTex printing flag.
        """
        for res_num, res in enumerate(results):
            var_sym = res[0]
            lfsr = res[3]
            cycle = res[1]
            initials = res[2]
            var_gen = lambdify((), var_sym, modules="mpmath")
            a_ = create_series_from_shift_reg(lfsr, initials, self.depth)
            b_ = (cycle * (self.depth // len(cycle)))[:self.depth]
            gcf = GeneralizedContinuedFraction(a_, b_)
            rate = calculate_convergence(gcf, lambdify((), var_sym, 'mpmath')())
            if not latex:
                print(str(res_num))
                print('lhs: ')
                sympy.pprint(var_sym)
                print('rhs :')
                gcf.print(8)
                print('lhs value: ' + mpmath.nstr(var_gen(), 50))
                print('rhs value: ' + mpmath.nstr(gcf.evaluate(), 50))
                print('a_n LFSR: {},\n With initialization: {}'.format(lfsr, initials))
                print('b_n period: ' + str(cycle))
                print("Converged with a rate of {} digits per term".format(mpmath.nstr(rate, 5)))
            else:
                equation = sympy.Eq(var_sym, gcf.sym_expression(5))
                print('\n\n' + str(res_num + 1) + '. $$ ' + sympy.latex(equation) + ' $$\n')
                print('$\\{a_n\\}$ LFSR: \\quad \\quad \\quad \\quad \\;' + str(lfsr))
                print('\n$\\{a_n\\}$ initialization: \\quad \\; ' + str(initials))
                print('\n$\\{b_n\\}$ Sequence period: \\! ' + str(cycle))
                print('\nConvergence rate: ${}$ digits per term\n\n'.format(mpmath.nstr(rate, 5)))

    def find_hits(self):
        """
        Use search engine to find results.
        :return: List of verified results, alongside a dictionary of similar results (of same numeric value).
        (The duplicates might prove useful later if we can find different sign series leading to different a series for
        same variation.)
        """
        with mpmath.workdps(self.enum_dps):
            start = time()
            # Search
            results = self.find_signed_rcf_conj()
            end = time()
            if self.do_print:
                print('That took {}s'.format(end - start))
        with mpmath.workdps(self.verify_dps):
            if self.do_print:
                print('Starting to verify results...')
            start = time()
            # Validate
            verified_results, recurring_value_results = self.verify_results(results)
            end = time()
            if self.do_print:
                print('{} results were verified.\nThat took {}'.format(len(verified_results), end - start))
            # Print if requested:
            if self.do_print:
                self.print_results(verified_results)
        return verified_results, recurring_value_results


def esma_search_wrapper(constant, custom_enum, poly_deg, coeff_lim,
                   cycle_range, min_deg, depth, out_dir=None, do_print=True):
    """
    A Wrapper for searching using ESMA (currently without multiprocessing).
    :param constant: sympy constant
    :param custom_enum: A ready-made enumeration that only requires substituting a variable 'x' with the constant. (opt)
    :param poly_deg: Maximum degree of numerator and denominator polynomials in the rational LHS.
    :param coeff_lim: Range of coefficients for the rational function on the LHS.
    :param cycle_range: Range of lengths for the sign sequence's period.
    :param min_deg: Used to exclude lower degree numerator and denominator polynomials in rational LHS. (opt)
    :param depth: Number of elements of a series to extract. Relates to length of typical LFSRs  of the consant. (opt)
    :param out_dir: Path of director to save result binaries. (opt)
    :param do_print: Print outputs (Used as False primarily for unit tests). (opt)
    :return: A list of results of the form [lhs(sympy), sign_period, a_initialization, a_LFSR].
             Dictionary, maps strings of values to lists of recurring results sharing value. (result format as above).
    """
    if depth is not None:
        enum = SignedRcfEnumeration(sym_constant=constant, cycle_len_range=cycle_range, depth=depth,
                                    coefficients_limit=coeff_lim, poly_deg=poly_deg, min_deg=min_deg,
                                    custom_enum=custom_enum, do_print=do_print)
    else:
        enum = SignedRcfEnumeration(sym_constant=constant, cycle_len_range=cycle_range, coefficients_limit=coeff_lim,
                                    poly_deg=poly_deg, min_deg=min_deg, custom_enum=custom_enum, do_print=do_print)
    result_list, recurring_results_dict = enum.find_hits()
    if out_dir:
        path = out_dir
        if not os.path.isdir(path):
            os.makedirs(path)
        dup = '/'.join([path, 'recurring_by_value_0'])
        res = '/'.join([path, 'res_list_0'])
        i = 1
        while os.path.exists(res):
            res = res[:-1] + str(i)
            dup = dup[:-1] + str(i)
            i += 1
        with open(res, 'wb') as f:
            pickle.dump(result_list, f)
        with open(dup, 'wb') as f:
            pickle.dump(recurring_results_dict, f)
    return result_list, recurring_results_dict


if __name__ =='__main__':
    print("Keep it simple, use the API")
