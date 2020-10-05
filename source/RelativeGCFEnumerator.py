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

from AbstractGCFEnumerator import *

# intermediate result - coefficients of lhs transformation, and compact polynomials for seeding an and bn series.
Match = namedtuple('Match', 'lhs_key rhs_an_poly rhs_bn_poly')
RefinedMatch = namedtuple('Match', 'lhs_key rhs_an_poly rhs_bn_poly lhs_match_idx')
FormattedResult = namedtuple('FormattedResult', 'LHS RHS GCF')
IterationMetadata = namedtuple('IterationMetadata', 'an_coef bn_coef iter_counter')

def gcf_calculation_to_precision(an_iterator, bn_iterator, result_precision,
    max_iters = 400):
    """
    Calculates the GCF value until all results are equal when regarding the precision required.
    If this isn't possible, compute until hitting max_iters 
    """
    # To reduce the number of divisions executed, we will loop for BURST_NUMBER iterations for 
    # num/denum, and then calculate the GCF (divide them). This will reduce the number of divisions
    # that take a lot of time
    # TODO - exponential backoff?
    BURST_NUMBER = 15
    computed_values = []

    prev_q = 0
    q = 1
    prev_p = 1
    # This is a ugly hack but it works. a[0] is handled before the rest here:
    p = an_iterator.__next__() # will place a[0] to p
    if p == 0:
        raise ZeroInAn()

    for i, (a_i_1, b_i) in enumerate(zip(an_iterator, bn_iterator)):
        if a_i_1 == 0:
            raise ZeroInAn()

        # a_i_1 is the (i+1)'th item of an, and b_i the the i'th item of bn
        tmp_a = q
        tmp_b = p

        q = a_i_1 * q + b_i * prev_q
        p = a_i_1 * p + b_i * prev_p
        
        prev_q = tmp_a
        prev_p = tmp_b

        # This can be much more efficient if divided to loops
        # save on calculations later
        items_computed = 0 
        if i % BURST_NUMBER:
            if q != 0:  # safety check
                computed_values.append(mpmath.mpf(p) / mpmath.mpf(q))
                items_computed += 1
            else:
                computed_values.append(mpmath.mpf(0))
                items_computed += 1

            # checking if the value stabilized
            # after the first 10 iterations, start checking the the last two values, which are (perhaps) the max and min, are equal considering the precession required
            if items_computed >= 2: # not in the first iteration
                if int(computed_values[-1]*result_precision) == \
                    int(computed_values[-2]*result_precision):
                    if i>50:
                        print(f'WOW super slow! {i}')
                    return computed_values[-1]
        
    # GCF didn't converge int time. guessing the avg between the 
    # last two calculations    
    avg = (computed_values[-1] + computed_values[-2])/2

    return avg


class RelativeGCFEnumerator(AbstractGCFEnumerator):
    """
        This enumerator calculates the gcf to an arbitrary dept, until getting to
        a stable value withing the precession required bounds.
        Useful for GCF that converges slowly.
    """

    def __init__(self, *args, **kwargs):
        print('using relative enumerator')
        super().__init__(*args, **kwargs)

    def _gcf_series_cached_iters(self, poly_a_domain, poly_b_domain, max_iters=100):
        """
        Calculates the approximate size for each domain, and using a nested loop over both
        domains return two iterators for each combination of a and b polynomials, where
        the smaller one is cached using CachedPolySeriesCalculator. 
        Also returns metadata about the execution stage, according to print_interval_percentage. 
        if its not the time to print stats, will return none for this variable 
        """
        size_a = self.get_an_length(poly_a_domain)
        size_b = self.get_bn_length(poly_b_domain)

        a_coef_iter = self.get_an_iterator(poly_a_domain)
        b_coef_iter = self.get_bn_iterator(poly_b_domain)

        iter_counter = 0
        # # we want to loop through all available combinations
        # if we'll just iter through both items, each internal series will be calculated
        # again for each external loop. Caching the smaller one and using it as the internal loop
        print(f'size a - {size_a} | size b - {size_b}')

        # This distinction is meant to make sure that we cache the smaller group
        if size_a > size_b:  # cache {bn} in RAM, iterate over an
            print('caching bn')
            bn_family = CachedPolySeriesCalculator()

            # external loop is not cached, and is for an
            
            for a_coef in a_coef_iter:
                a_degree, a_leading_coef = get_poly_deg_and_leading_coef(a_coef)
                an_cache = CachedSeriesCalculator(a_coef)

                for b_coef, bn_iterator in bn_family.iter_family(b_coef_iter, max_iters):
                    iter_counter += 1
                    
                    yield \
                        an_cache.iter_series_items(max_iters), \
                        bn_iterator, \
                        IterationMetadata(a_coef, b_coef, iter_counter)

        else:  # cache {an} in RAM, iterate over bn
            print('caching an')
            an_family = CachedPolySeriesCalculator()

            # external loop is not cached, and is for an
            for b_coef in b_coef_iter:
                bn_cache = CachedSeriesCalculator(b_coef)
                b_degree, b_leading_coef = get_poly_deg_and_leading_coef(b_coef)

                for a_coef, an_iterator in an_family.iter_family(a_coef_iter, max_iters):
                    iter_counter += 1
                    yield \
                        an_iterator, \
                        bn_cache.iter_series_items(max_iters), \
                        IterationMetadata(a_coef, b_coef, iter_counter)


    def _first_enumeration(self, poly_a: List[List], poly_b: List[List], print_results: bool):
        size_a = self.get_an_length(poly_a)
        size_b = self.get_bn_length(poly_b)

        num_iterations = size_b * size_a

        start = time()
        key_factor = 1 / self.threshold

        results = []  # list of intermediate results        

        for an_iter, bn_iter, metadata in self._gcf_series_cached_iters(poly_a, poly_b, g_N_initial_search_terms):

            try:
                gcf_val = gcf_calculation_to_precision(an_iter, bn_iter, g_N_initial_key_length)
            except ZeroInAn:
                continue

            key = int(gcf_val * key_factor)
            if key in self.hash_table:  # find hits in hash table
                results.append(Match(key, metadata.an_coef, metadata.bn_coef))

            # This is a fast way to check if iter counter hit some print interval 
            if metadata.iter_counter & 0x10000 == metadata.iter_counter & 0x1ffff:  # print status.
                print(
                    f'passed {metadata.iter_counter} out of {num_iterations} ({round(100. * metadata.iter_counter / num_iterations, 2)}%). found so far {len(results)} results')
                print(f'currently at an = {metadata.an_coef} bn = {metadata.bn_coef}')
 
        if print_results:
            print(f'created results after {time() - start}s')
            print(f'found {len(results)} results')
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