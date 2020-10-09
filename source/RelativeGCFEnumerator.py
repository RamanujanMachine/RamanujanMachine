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
    min_iters=1, burst_number=15):
    """
    Calculates the GCF value until all results are equal when regarding the precision required.
    If this isn't possible, compute until hitting max_iters 
    """
    # To reduce the number of divisions executed, we will loop for BURST_NUMBER iterations for 
    # num/denum, and then calculate the GCF (divide them). This will reduce the number of divisions
    # that take a lot of time
    # TODO - exponential backoff?
    computed_values = []
    items_computed = 0 

    # Some GCFs converge by oscillation. its better to take odd burst_number to test
    # values from both sub-series
    burst_number = burst_number if burst_number % 2 ==1 else burst_number +1
    next_gcf_calculation = burst_number if burst_number >= min_iters else min_iters
    
    precision_factor = 10**result_precision
    prev_q = 0
    q = 1
    prev_p = 1
    # This is a ugly hack but it works. a[0] is handled before the rest here:
    p = an_iterator.__next__() # will place a[0] to p
    bn_iterator.__next__() # b0 is discarded

    if p == 0:
        raise ZeroInAn()

    for i, (a_i, b_i) in enumerate(zip(an_iterator, bn_iterator)):
        if a_i == 0:
            raise ZeroInAn()

        # a_i is the (i+1)'th item of an, and b_i the the i'th item of bn
        tmp_a = q
        tmp_b = p

        q = a_i * q + b_i * prev_q
        p = a_i * p + b_i * prev_p
        
        prev_q = tmp_a
        prev_p = tmp_b

        # using module so many times take noticeable resorces. using this method 
        # to decide when to divide next should save some time
        if i == next_gcf_calculation:
            next_gcf_calculation += burst_number
            if q != 0:  # safety check
                computed_values.append(mpmath.mpf(p) / mpmath.mpf(q))
                items_computed += 1
            else:
                computed_values.append(mpmath.mpf(0))
                items_computed += 1

            # checking if the value stabilized
            # after the first 10 iterations, start checking the the last two values, which are (perhaps) the max and min, are equal considering the precession required
            if items_computed >= 2: # not in the first iteration
                if int(computed_values[-1]*precision_factor) == \
                    int(computed_values[-2]*precision_factor):
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

        Note-
            Currently using costume series generators is not allowed, since those 
            operate by computing a fixed n items, as opposed to the fluid approach
            which suggested here. To implement it, one need to make those an iterative
            function, and pass it to cached_series_iterator (and also add support there)
    """

    def __init__(self, *args, **kwargs):
        print('using relative enumerator')
        super().__init__(*args, **kwargs)

    def _gcf_series_cached_iters(self, poly_domains, max_iters=100):
        """
        Calculates the approximate size for each domain, and using a nested loop over both
        domains return two iterators for each combination of a and b polynomials, where
        the smaller one is cached using CachedPolySeriesCalculator. 
        Also returns metadata about the execution stage, according to print_interval_percentage. 
        if its not the time to print stats, will return none for this variable 
        """
        size_a = poly_domains.an_length
        size_b = poly_domains.bn_length

        a_coef_iter, b_coef_iter = poly_domains.get_individual_polys_generators()

        an_series_iter, bn_series_iter = poly_domains.get_calculation_method()

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
                #a_degree, a_leading_coef = get_poly_deg_and_leading_coef(a_coef)
                #an_cache = CachedSeriesCalculator(a_coef)

                b_coef_iter = poly_domains.get_b_coef_iterator()
                for b_coef, bn_iterator in bn_family.iter_family(b_coef_iter, max_iters,
                    series_iterator=bn_series_iter):
                    iter_counter += 1

                    yield \
                        an_series_iter(a_coef, max_iters, 0), \
                        bn_iterator, \
                        IterationMetadata(a_coef, b_coef, iter_counter)

        else:  # cache {an} in RAM, iterate over bn
            print('caching an')
            an_family = CachedPolySeriesCalculator()

            # external loop is not cached, and is for an
            for b_coef in b_coef_iter:
                #bn_cache = CachedSeriesCalculator(b_coef)
                #b_degree, b_leading_coef = get_poly_deg_and_leading_coef(b_coef)
                
                a_coef_iter = poly_domains.get_a_coef_iterator()
                for a_coef, an_iterator in an_family.iter_family(a_coef_iter, max_iters,
                    series_iterator=an_series_iter):
                    iter_counter += 1
                    yield \
                        an_iterator, \
                        bn_series_iter(b_coef, max_iters, 0), \
                        IterationMetadata(a_coef, b_coef, iter_counter)


    def _first_enumeration(self, poly_domains, print_results: bool):
        num_iterations = poly_domains.num_iterations

        start = time()
        key_factor = 1 / self.threshold

        results = []  # list of intermediate results        

        for an_iter, bn_iter, metadata in self._gcf_series_cached_iters(poly_domains, 1000):
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
                    f'passed {metadata.iter_counter} out of {poly_domains.num_iterations} ' + \
                    f'({round(100. * metadata.iter_counter / poly_domains.num_iterations, 2)}%). ' + \
                    f' found so far {len(results)} results')
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

            # # create a_n, b_n with huge length, calculate gcf, and verify result.
            # an = self.create_an_series(res.rhs_an_poly, g_N_verify_terms)
            # bn = self.create_bn_series(res.rhs_bn_poly, g_N_verify_terms)
            # gcf = EfficientGCF(an, bn)
            an_iter_func, bn_iter_func = self.poly_domains_generator.get_calculation_method()
            an_iter = an_iter_func(res.rhs_an_poly, 100000, start_n=0)
            bn_iter = bn_iter_func(res.rhs_bn_poly, 100000, start_n=0)

            gcf = gcf_calculation_to_precision(an_iter, bn_iter, g_N_verify_compare_length, 
                min_iters=2000, burst_number=500)
            rhs_str = mpmath.nstr(gcf, g_N_verify_compare_length)

            for i, val in enumerate(all_matches):
                val_str = mpmath.nstr(val, g_N_verify_compare_length)
                if val_str == rhs_str:
                    # This patch is ment to allow support for multiple matches for an
                    # LHS key, i will later be used to determind which item in the LHS dict
                    # was matched
                    results.append(RefinedMatch(*res, i))

        return results