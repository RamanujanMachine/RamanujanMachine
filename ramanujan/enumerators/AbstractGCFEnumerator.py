from time import time
from typing import List
from collections import namedtuple
from collections.abc import Iterable
from sympy import lambdify
from abc import ABCMeta, abstractmethod

from ramanujan.utils.mobius import GeneralizedContinuedFraction
from ramanujan.utils.utils import find_polynomial_series_coefficients, create_mpf_const_generator, \
    get_series_items_from_iter
from ramanujan.utils.convergence_rate import calculate_convergence
from ramanujan.constants import *

Match = namedtuple('Match', 'lhs_key rhs_an_poly rhs_bn_poly')
RefinedMatch = namedtuple('RefinedMatch', 'lhs_key rhs_an_poly rhs_bn_poly lhs_match_idx c_top c_bot')
FormattedResult = namedtuple('FormattedResult', 'LHS RHS GCF')


def get_size_of_nested_list(list_of_elem):
    """ Get number of elements in a nested list"""
    count = 0
    # Iterate over the list
    for elem in list_of_elem:
        # Check if type of element is list
        if isinstance(elem, Iterable):
            # Again call this function to get the size of this element
            count += get_size_of_nested_list(elem)
        else:
            count += 1
    return count


class ZeroInAn(Exception):
    pass


class AbstractGCFEnumerator(metaclass=ABCMeta):
    """
        This is an abstract class for RHS searches 'engines'.

        basically, this is a 2 step procedure:
        1) first enumeration - enumerate over all rhs combinations, find hits in lhs hash table.
            For each hit, calculate the GCF to a higher dept.
        2) refine results - take results from (1) and validate them to 100 decimal digits.

        Functions implemented in this abstract class will translate the given constants to
        execution limits (e.g. cast search limit to convergence threshold), handle printing
        results, and making mpmath work to a sufficient precision on every location

        Enumerators should implement the following:
            find_initial_hits
            refine_results
    """

    def __init__(self, hash_table, poly_domains, sym_constants):
        """
        initialize search engine.
        :param hash_table: LHSHashTable object storing the constant's permutations. Used for
            querying computed values for
        :param poly_domains: An poly_domain object that will generate polynomials to iter through, and
            supply functions for calculating items in each polynomial given
        :param sym_constants: sympy constants
        """
        # constants
        self.threshold = 1 * 10 ** (-g_N_initial_key_length)  # key length
        self.enum_dps = g_N_initial_search_dps  # working decimal precision for first enumeration
        self.verify_dps = g_N_verify_dps  # working decimal precision for validating results.
        self.const_sym = sym_constants
        self.constants_generator = create_mpf_const_generator(sym_constants)

        # expand poly domains object
        # there are two methods to generate and iter over domains.  the newer one uses poly_domains only,
        # but the old one still uses the rest of the arguments.
        # generating them here to avoid breaking older enumerators
        self.poly_domains = poly_domains
        a_iterator_func, b_iterator_func = poly_domains.get_calculation_method()
        self.create_an_series = \
            lambda coefs, items: get_series_items_from_iter(a_iterator_func, coefs, items)
        self.create_bn_series = \
            lambda coefs, items: get_series_items_from_iter(b_iterator_func, coefs, items)
        self.get_an_length = poly_domains.get_an_length
        self.get_bn_length = poly_domains.get_bn_length

        self.get_an_iterator = poly_domains.get_a_coef_iterator
        self.get_bn_iterator = poly_domains.get_b_coef_iterator

        # store lhs_hash_table
        self.hash_table = hash_table

    def __get_formatted_results(self, results: List[RefinedMatch]) -> List[FormattedResult]:
        ret = []
        for r in results:
            an = self.create_an_series(r.rhs_an_poly, 250)
            bn = self.create_bn_series(r.rhs_bn_poly, 250)
            print_length = max(max(get_size_of_nested_list(r.rhs_an_poly), get_size_of_nested_list(r.rhs_bn_poly)), 5)
            gcf = GeneralizedContinuedFraction(an, bn[1:])
            sym_lhs = self.hash_table.evaluate_sym(r.lhs_key, self.const_sym)[r.lhs_match_idx]
            ret.append(FormattedResult(sym_lhs, gcf.sym_expression(print_length), gcf))
        return ret

    def __get_formatted_polynomials(self, result: RefinedMatch):
        def sym_poly(poly_deg, poly_terms):
            poly = list(reversed(find_polynomial_series_coefficients(poly_deg, poly_terms, 0)))
            n = sympy.Symbol('n')
            poly_sym = 0
            for i in range(len(poly)):
                poly_sym += n ** i * poly[i]
            return poly_sym

        an_poly_max_deg = get_size_of_nested_list(result.rhs_an_poly)
        an = self.create_an_series(result.rhs_an_poly, an_poly_max_deg + 1)
        bn_poly_max_deg = get_size_of_nested_list(result.rhs_bn_poly)
        bn = self.create_bn_series(result.rhs_bn_poly, bn_poly_max_deg + 1)
        an_eq = sympy.Eq(sympy.Symbol('a(n)'), sym_poly(an_poly_max_deg, an))
        bn_eq = sympy.Eq(sympy.Symbol('b(n)'), sym_poly(bn_poly_max_deg, bn))
        return an_eq, bn_eq

    def print_results(self, results: List[RefinedMatch], formatting='unicode', convergence_rate=True):
        """
        pretty print the the results.
        :param convergence_rate: if True calculate convergence rate and print it as well.
        :param results: list of final results as received from refine_results.
        :param formatting: allowed print formats are 'unicode' and 'latex'
        """
        allowed_formats = ['unicode', 'latex']
        formatted_results = self.__get_formatted_results(results)
        if formatting not in allowed_formats:
            print("unknown format, allowed formats are: {}".format(allowed_formats))
            return
        for r, raw_r in zip(formatted_results, results):
            result = sympy.Eq(r.LHS, r.RHS)
            if formatting == 'latex':
                print(f'$$ {sympy.latex(result)} $$')
                print(f'$$ {sympy.latex(self.__get_formatted_polynomials(raw_r))} $$\n')
            else:
                sympy.pprint(result)
                sympy.pprint(self.__get_formatted_polynomials(raw_r))
                print('')
            if convergence_rate:
                with mpmath.workdps(self.verify_dps):
                    rate = calculate_convergence(r.GCF, lambdify((), r.LHS, 'mpmath')())
                print("Converged with a rate of {} digits per term".format(mpmath.nstr(rate, 5)))

    def convert_results_to_latex(self, results: List[RefinedMatch]):
        results_in_latex = []
        formatted_results = self.__get_formatted_results(results)
        for r in formatted_results:
            equation = sympy.Eq(r.LHS, r.RHS)
            results_in_latex.append(sympy.latex(equation))
        return results_in_latex

    def find_initial_hits(self, verbose=True):
        """
        use search engine to find results (steps (1) explained in __init__ docstring)
        :param verbose: if true, pretty print results at the end.
        :return: initial results results.
        """
        with mpmath.workdps(self.enum_dps):
            if verbose:
                print('starting preliminary search...')
            start = time()
            # step (2)
            intermediate_results = self._first_enumeration(verbose)
            end = time()
            if verbose:
                print(f'that took {end - start}s')

        with mpmath.workdps(self.verify_dps * 2):
            if verbose:
                print('calculating intermediate results to a higher precision...')
            start = time()
            results = self._improve_results_precision(intermediate_results, verbose)
            end = time()
            if verbose:
                print(f'that took {end - start}s')

        return results

    @abstractmethod
    def _first_enumeration(self, verbose: bool):
        # override by child
        pass

    def _improve_results_precision(self, intermediate_results, verbose: bool):
        """
        Calculates intermediate results GCFs to a higher dept, yielding more precise results.
        """
        pass

    def refine_results(self, results):
        with mpmath.workdps(self.verify_dps * 2):
            print('starting to verify results...')
            start = time()
            refined_results = self._refine_results(results, True)  # step (3)
            end = time()
            print(f'that took {end - start}s')
        return refined_results

    @abstractmethod
    def _refine_results(self, intermediate_results: List[Match], verbose=True):
        # override by child
        pass

    def full_execution(self):
        first_iteration = self.find_initial_hits()
        refined_results = self.refine_results(first_iteration)
        return refined_results
