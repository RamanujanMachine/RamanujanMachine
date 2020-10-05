import os
import pickle
#import itertools
#import multiprocessing
#import struct
#from functools import partial, reduce
#from datetime import datetime
from time import time
#from math import gcd
from typing import List, Iterator, Callable
from collections import namedtuple
from collections.abc import Iterable
import mpmath
import sympy
from sympy import lambdify
from latex import generate_latex
from pybloom_live import BloomFilter
from mobius import GeneralizedContinuedFraction, EfficientGCF
from convergence_rate import calculate_convergence
from series_generators import SeriesGeneratorClass, CartesianProductAnGenerator, CartesianProductBnGenerator
from utils import find_polynomial_series_coefficients, get_poly_deg_and_leading_coef
from LHSHashTable import LHSHashTable
from cached_poly_series_calculator import CachedPolySeriesCalculator
from cached_series_calculator import CachedSeriesCalculator
from series_generators import iter_series_items_from_compact_poly
from abc import ABCMeta, abstractmethod

# intermediate result - coefficients of lhs transformation, and compact polynomials for seeding an and bn series.
Match = namedtuple('Match', 'lhs_key rhs_an_poly rhs_bn_poly')
RefinedMatch = namedtuple('Match', 'lhs_key rhs_an_poly rhs_bn_poly lhs_match_idx')
FormattedResult = namedtuple('FormattedResult', 'LHS RHS GCF')
IterationMetadata = namedtuple('IterationMetadata', 'an_coef bn_coef iter_counter')

# global hyper-parameters
g_N_verify_terms = 1000  # number of CF terms to calculate in __refine_results. (verify hits)
g_N_verify_compare_length = 100  # number of digits to compare in __refine_results. (verify hits)
g_N_verify_dps = 2000  # working decimal precision in __refine_results. (verify hits)
g_N_initial_search_terms = 32  # number of CF terms to calculate in __first_enumeration (initial search)
g_N_initial_key_length = 10  # number of digits to compare in __first_enumeration (initial search)
g_N_initial_search_dps = 50  # working decimal precision in __refine_results. (verify hits)


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

        basically, this is a 3 step procedure:
        1) load / initialize lhs hash table.
        2) first enumeration - enumerate over all rhs combinations, find hits in lhs hash table.
        3) refine results - take results from (2) and validate them to 100 decimal digits.
        
        Functions in this class handles initiating the class with
        constants, the hash table, and printing results
        
        Enumerators should implement the following:
            find_initial_hits
            refine_results

    """
    def __init__(self, sym_constants, lhs_search_limit, saved_hash,
                 an_generator: SeriesGeneratorClass = CartesianProductAnGenerator(),
                 bn_generator: SeriesGeneratorClass = CartesianProductBnGenerator()):
        """
        initialize search engine.
        :param sym_constants: sympy constants
        :param lhs_search_limit: range of coefficients for left hand side.
        :param saved_hash: path to saved hash.
        :param an_generator: generating function for {an} series
        :param bn_generator: generating function for {bn} series
        """
        self.threshold = 1 * 10 ** (-g_N_initial_key_length)  # key length
        self.enum_dps = g_N_initial_search_dps  # working decimal precision for first enumeration
        self.verify_dps = g_N_verify_dps  # working decimal precision for validating results.
        self.lhs_limit = lhs_search_limit
        self.const_sym = sym_constants
        self.constants_generator = []
        for i in range(len(sym_constants)):
            try:
                self.constants_generator.append(lambdify((), sym_constants[i], modules="mpmath"))
            except AttributeError:  # Hackish constant
                self.constants_generator.append(sym_constants[i].mpf_val)

        self.create_an_series = an_generator.get_function()
        self.get_an_length = an_generator.get_num_iterations
        self.get_an_iterator = an_generator.get_iterator
        self.create_bn_series = bn_generator.get_function()
        self.get_bn_length = bn_generator.get_num_iterations
        self.get_bn_iterator = bn_generator.get_iterator

        if not os.path.isfile(saved_hash):
            print('no previous hash table given, initializing hash table...')
            with mpmath.workdps(self.enum_dps):
                constants = [const() for const in self.constants_generator]
                start = time()
                self.hash_table = LHSHashTable(
                    saved_hash,
                    self.lhs_limit,
                    constants,  # constant
                    self.threshold)  # length of key
                end = time()
                print(f'that took {end - start}s')
        else:
            self.hash_table = LHSHashTable.load_from(saved_hash)

    def __get_formatted_results(self, results: List[Match]) -> List[FormattedResult]:
        ret = []
        for r in results:
            an = self.create_an_series(r.rhs_an_poly, 250)
            bn = self.create_bn_series(r.rhs_bn_poly, 250)
            print_length = max(max(get_size_of_nested_list(r.rhs_an_poly), get_size_of_nested_list(r.rhs_bn_poly)), 5)
            gcf = GeneralizedContinuedFraction(an, bn)
            sym_lhs = self.hash_table.evaluate_sym(r.lhs_key, self.const_sym)[r.lhs_match_idx]
            ret.append(FormattedResult(sym_lhs, gcf.sym_expression(print_length), gcf))
        return ret

    def __get_formatted_polynomials(self, result: Match):
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

    def print_results(self, results: List[Match], latex=False, convergence_rate=True):
        """
        pretty print the the results.
        :param convergence_rate: if True calculate convergence rate and print it as well.
        :param results: list of final results as received from refine_results.
        :param latex: if True print in latex form, otherwise pretty print in unicode.
        """
        formatted_results = self.__get_formatted_results(results)
        for r, raw_r in zip(formatted_results, results):
            result = sympy.Eq(r.LHS, r.RHS)
            if latex:
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

    def convert_results_to_latex(self, results: List[Match]):
        results_in_latex = []
        formatted_results = self.__get_formatted_results(results)
        for r in formatted_results:
            equation = sympy.Eq(r.LHS, r.RHS)
            results_in_latex.append(sympy.latex(equation))
        return results_in_latex

    def find_initial_hits(self, poly_a: List[List], poly_b: List[List], print_results=True):
        """
        use search engine to find results (steps (2) and (3) explained in __init__ docstring)
        :param poly_a: explained in docstring of __first_enumeration
        :param poly_b: explained in docstring of __first_enumeration
        :param print_results: if true, pretty print results at the end.
        :return: initial results results.
        """
        with mpmath.workdps(self.enum_dps):
            if print_results:
                print('starting preliminary search...')
            start = time()
            # step (2)
            results = self._first_enumeration(poly_a, poly_b, print_results)
            end = time()
            if print_results:
                print(f'that took {end - start}s')
        return results

    @abstractmethod    
    def _first_enumeration(self, poly_a: List[List], poly_b: List[List], print_results: bool):
        # override by child
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
    def _refine_results(self, intermediate_results: List[Match], print_results=True):
        # override by child
        pass 
    