import os
import pickle
import itertools
import multiprocessing
import struct
from functools import partial, reduce
from datetime import datetime
from time import time
from math import gcd
from typing import List, Iterator, Callable
from collections import namedtuple
from collections.abc import Iterable
import mpmath
import sympy
from sympy import lambdify
from pybloom_live import BloomFilter

from ramanujan.utils.latex import generate_latex
from ramanujan.mobius import GeneralizedContinuedFraction, EfficientGCF
from ramanujan.utils.convergence_rate import calculate_convergence
from ramanujan.series_generators import SeriesGeneratorClass, CartesianProductAnGenerator, CartesianProductBnGenerator
from ramanujan.utils.utils import find_polynomial_series_coefficients
from ramanujan.LHSHashTable import LHSHashTable

from .AbstractGCFEnumerator import *

class EfficentGCFEnumerator(AbstractGCFEnumerator):
    """
        This enumerator maximizes on efficiency for calculating GCF to a
        precise dept. will produce false results when calculating GCFs
        that converges slowly
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def __create_series_list(coefficient_iter: Iterator,
                             series_generator: Callable[[List[int], int], List[int]],
                             filter_from_1=False) -> [List[int], List[int]]:
        coef_list = list(coefficient_iter)
        # create a_n and b_n series fro coefficients.
        series_list = [series_generator(coef_list[i], g_N_initial_search_terms) for i in range(len(coef_list))]
        # filter out all options resulting in '0' in any series term.
        if filter_from_1:
            series_filter = [0 not in an[1:] for an in series_list]
        else:
            series_filter = [0 not in an for an in series_list]
        series_list = list(itertools.compress(series_list, series_filter))
        coef_list = list(itertools.compress(coef_list, series_filter))
        return coef_list, series_list

    def _first_enumeration(self, poly_a: List[List], poly_b: List[List], print_results: bool):
        """
        this is usually the bottleneck of the search.
        we calculate general continued fractions of type K(bn,an). 'an' and 'bn' are polynomial series.
        these polynomials take the form of n(n(..(n*c_1 + c_0) + c_2)..)+c_k.
        poly parameters are a list of coefficients c_i. then the enumeration takes place on all possible products.
        for example, if poly_a is [[1,2],[2,3]], then the product polynomials are:
           possible [c0,c1] = { [1,2] , [1,3], [2,2], [2,3] }.
        this explodes exponentially -
        example: fs poly_a.shape = poly_b.shape = [3,5]     (2 polynomials of degree 2), then the number of
        total permutations is: (a poly possibilities) X (b poly possibilities) = (5X5X5) X (5X5X5) = 5**6
        we search on all possible gcf with polynomials defined by parameters, and try to find hits in hash table.
        :param poly_a: compact polynomial form of 'an' (list of lists)
        :param poly_b: compact polynomial form of 'an' (list of lists)
        :param print_results: if True print the status of calculation.
        :return: intermediate results (list of 'Match')
        """

        def efficient_gcf_calculation():
            """
            enclosure. a_, b_, and key_factor are used from outer scope.
            moved here from mobius.EfficientGCF to optimize performance.
            :return: key for LHS hash table
            """
            prev_q = 0
            q = 1
            prev_p = 1
            p = a_[0]
            for i in range(1, len(a_)):
                tmp_a = q
                tmp_b = p
                q = a_[i] * q + b_[i - 1] * prev_q
                p = a_[i] * p + b_[i - 1] * prev_p
                prev_q = tmp_a
                prev_p = tmp_b
            if q == 0:  # safety check
                value = 0
            else:
                value = mpmath.mpf(p) / mpmath.mpf(q)
            return int(value * key_factor)  # calculate hash key of gcf value

        start = time()
        a_coef_iter = self.get_an_iterator(poly_a)  # all coefficients possibilities for 'a_n'
        b_coef_iter = self.get_bn_iterator(poly_b)
        size_b = self.get_bn_length(poly_b)
        size_a = self.get_an_length(poly_a)
        num_iterations = size_b * size_a
        key_factor = 1 / self.threshold

        counter = 0  # number of permutations passed
        print_counter = counter
        results = []  # list of intermediate results

        if size_a > size_b:  # cache {bn} in RAM, iterate over an
            b_coef_list, bn_list = self.__create_series_list(b_coef_iter, self.create_bn_series)
            real_bn_size = len(bn_list)
            num_iterations = (num_iterations // self.get_bn_length(poly_b)) * real_bn_size
            if print_results:
                print(f'created final enumerations filters after {time() - start}s')
            start = time()
            for a_coef in a_coef_iter:
                an = self.create_an_series(a_coef, g_N_initial_search_terms)
                if 0 in an[1:]:  # a_0 is allowed to be 0.
                    counter += real_bn_size
                    print_counter += real_bn_size
                    continue
                for bn_coef in zip(bn_list, b_coef_list):
                    # evaluation of GCF: taken from mobius.EfficientGCF and moved here to avoid function call overhead.
                    a_ = an
                    b_ = bn_coef[0]
                    key = efficient_gcf_calculation()  # calculate hash key of gcf value

                    if key in self.hash_table:  # find hits in hash table
                        results.append(Match(key, a_coef, bn_coef[1]))
                    if print_results:
                        counter += 1
                        print_counter += 1
                        if print_counter >= 100000:  # print status.
                            print_counter = 0
                            print(
                                f'passed {counter} out of {num_iterations} ({round(100. * counter / num_iterations, 2)}%). found so far {len(results)} results')

        else:  # cache {an} in RAM, iterate over bn
            a_coef_list, an_list = self.__create_series_list(a_coef_iter, self.create_an_series, filter_from_1=True)
            real_an_size = len(an_list)
            num_iterations = (num_iterations // self.get_an_length(poly_a)) * real_an_size
            if print_results:
                print(f'created final enumerations filters after {time() - start}s')
            start = time()
            for b_coef in b_coef_iter:
                bn = self.create_bn_series(b_coef, g_N_initial_search_terms)
                if 0 in bn:
                    counter += real_an_size
                    print_counter += real_an_size
                    continue
                for an_coef in zip(an_list, a_coef_list):
                    a_ = an_coef[0]
                    b_ = bn
                    key = efficient_gcf_calculation()  # calculate hash key of gcf value

                    if key in self.hash_table:  # find hits in hash table
                        results.append(Match(key, an_coef[1], b_coef))
                    if print_results:
                        counter += 1
                        print_counter += 1
                        if print_counter >= 100000:  # print status.
                            print_counter = 0
                            print(
                                f'passed {counter} out of {num_iterations} ({round(100. * counter / num_iterations, 2)}%). found so far {len(results)} results')

        if print_results:
            print(f'created results after {time() - start}s')
        return results

    def _refine_results(self, intermediate_results: List[Match], print_results=True):
        """
        validate intermediate results to 100 digit precision
        :param intermediate_results:  list of results from first enumeration
        :param print_results: if true print status.
        :return: final results.
        """
        results = []
        counter = 0
        n_iterations = len(intermediate_results)
        constant_vals = [const() for const in self.constants_generator]
        for res in intermediate_results:
            counter += 1
            if (counter % 50) == 0 and print_results:
                print('passed {} permutations out of {}. found so far {} matches'.format(
                    counter, n_iterations, len(results)))
            try:
                all_matches = self.hash_table.evaluate(res.lhs_key, constant_vals)
                # check if all values enounter are not inf or nan
                # TODO - consider mpmath.isnormal(val)
                if not all([ not (mpmath.isinf(val) or mpmath.isnan(val)) for val in all_matches]):  # safety
                    print('Something wicked happend!')
                    print(f'Encountered a NAN or inf in LHS db, at {res.lhs_key}, {constant_vals}')
                    continue
            except (ZeroDivisionError, KeyError) as e:
                continue

            # create a_n, b_n with huge length, calculate gcf, and verify result.
            an = self.create_an_series(res.rhs_an_poly, g_N_verify_terms)
            bn = self.create_bn_series(res.rhs_bn_poly, g_N_verify_terms)
            gcf = EfficientGCF(an, bn)
            rhs_str = mpmath.nstr(gcf.evaluate(), g_N_verify_compare_length)
            
            for i, val in enumerate(all_matches):
                val_str = mpmath.nstr(val, g_N_verify_compare_length)
                if val_str == rhs_str:
                    # This patch is ment to allow support for multiple matches for an
                    # LHS key, i will later be used to determind which item in the LHS dict
                    # was matched
                    results.append(RefinedMatch(*res, i))

        return results