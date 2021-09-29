from time import time
from typing import List
from collections import namedtuple
import mpmath

from ramanujan.CachedSeries import CachedSeries
from ramanujan.constants import g_N_verify_compare_length, g_N_initial_key_length
from .AbstractGCFEnumerator import AbstractGCFEnumerator, Match
from ramanujan.utils.utils import trunc_division

RefinedMatch = namedtuple('RefinedMatch', 'lhs_key rhs_an_poly rhs_bn_poly lhs_match_idx c_top c_bot precision')

IterationMetadata = namedtuple('IterationMetadata', 'an_coef bn_coef')

FIRST_STEP_MIN_ITERS = 7
FIRST_STEP_MAX_ITERS = 100
FIRST_STEP_BURST_NUMBER = 7

SECOND_STEP_MIN_ITERS = 7
SECOND_STEP_MAX_ITERS = 10_000
SECOND_STEP_BURST_NUMBER = 7


class ZeroInAn(Exception):
    pass


class NotConverging(Exception):
    pass


def gcf_calculation_to_precision(an_iterator, bn_iterator, result_precision, min_iters, burst_number):
    """
    Calculate the GCF's value to a changing depth that depends on the GCF's convergence rate.

    Most GCFs will converge by oscillation. If we calculate the GCF's value on the n and n+1 depth, we can conclude
    (for most GCFs) that the resulting value will be between those two values. Using this principle, we'll stop the
    GCF's calculation when those two values have the same decimal value under the required approximation.

    If the value didn't converge and the iterators are not yielding any new items, we'll take the last two values
    calculated for the GCF and compare their digits. If the digits match - we'll assume that it is a real value.

    Since huge int division is a costly process, we'll do so only every burst_number of calculations.

    Returns an integer with the matching key for the GCF, which is int(gcf_value*precision_factor) and the number of
    digits calculated.
    """
    computed_values = []
    items_computed = 0

    # Some GCFs converge by oscillation. its better to take odd burst_number to test
    # values from both sub-series
    burst_number = burst_number if burst_number % 2 == 1 else burst_number + 1
    next_gcf_calculation = burst_number if burst_number >= min_iters else min_iters

    precision_factor = 10 ** result_precision
    prev_q = 0
    q = 1
    prev_p = 1
    
    p = next(an_iterator)  # will place a[0] to p
    next(bn_iterator)  # b0 is discarded

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

        if i == next_gcf_calculation:
            next_gcf_calculation += burst_number * (len(computed_values) + 1)
            if q != 0:  # safety check

                computed_values.append(trunc_division(precision_factor * p, q))
                items_computed += 1
            else:
                raise ZeroDivisionError()

            if items_computed >= 2:
                if computed_values[-1] == computed_values[-2]:
                    return computed_values[-1], result_precision

                if items_computed >= 3:
                    # if the GCF oscillates between two sub-series - we expect the distance between the two to 
                    # monotonically decrease. If the GCF doesn't oscillate, we expect the distance between two following
                    # values to monotonically decrease as well. If none of those apply, we halt the calculation here.
                    if abs(computed_values[-2] - computed_values[-1]) > abs(computed_values[-3] - computed_values[-2]):
                        raise NotConverging("Not converging")

    # we'll take the last two calculations, and check for matching digits.
    # Once the two doesn't match, we'll know that we cannot trust the following digits.
    res = ''
    for i, (c1, c2) in enumerate(zip(str(computed_values[-2]), str(computed_values[-1]))):
        if c1 != c2:
            break
        res += c1

    res += '0' * (len(str(computed_values[-1]))-i)

    return int(res), i
        

class RelativeGCFEnumerator(AbstractGCFEnumerator):
    """
        This enumerator calculates the GCF to an arbitrary depth, until getting to a stable value within the precision
        required bounds.
        Useful for GCF that converges slowly.
    """

    def __init__(self, *args, **kwargs):
        print('using relative enumerator')
        super().__init__(*args, **kwargs)

    def _iter_domains_with_cache(self, max_iters):
        """
        Yields an and bn pairs in the given domain while keeping cache for one of the series (an or bn)
        to avoid double calculations.

        We allow calculations up to max_iters, and we'll avoid calculating items we don't use, so the 
        item in an or bn is only calculated when it's first yielded and then cached.
        """
        size_a = self.poly_domains.an_length
        size_b = self.poly_domains.bn_length

        an_series_iter, bn_series_iter = self.poly_domains.get_calculation_method()

        # The series on the outer loop is only used on one iteration of the out loop. So
        # we'll cache only the current series for it. All of the domains for the inner series are cached
        if size_a > size_b:  # cache bn
            bn_cache = {}
            an_cache = CachedSeries((0,), an_series_iter)
            for an_coefs, bn_coefs in self.poly_domains.iter_polys(primary_looped_domain='a'):
                if bn_coefs not in bn_cache:
                    bn_cache[bn_coefs] = CachedSeries(bn_coefs, bn_series_iter)
                if an_cache.poly_coefs != an_coefs:
                    an_cache = CachedSeries(an_coefs, an_series_iter)

                yield (
                    an_cache.iter_series_items(max_iters=max_iters),
                    bn_cache[bn_coefs].iter_series_items(max_iters=max_iters),
                    IterationMetadata(an_coefs, bn_coefs))

        else:  # cache an
            an_cache = {}
            bn_cache = CachedSeries((0,), bn_series_iter)
            for an_coefs, bn_coefs in self.poly_domains.iter_polys(primary_looped_domain='b'):
                if an_coefs not in an_cache:
                    an_cache[an_coefs] = CachedSeries(an_coefs, an_series_iter)
                if bn_cache.poly_coefs != bn_coefs:
                    bn_cache = CachedSeries(bn_coefs, bn_series_iter)

                yield (
                    an_cache[an_coefs].iter_series_items(max_iters=max_iters),
                    bn_cache.iter_series_items(max_iters=max_iters),
                    IterationMetadata(an_coefs, bn_coefs))

    def _first_enumeration(self, verbose: bool):
        """
        Calculate the GCD to a low precision and check for hits with the bloom filter.
        """
        start = time()

        results = []  # list of intermediate results        
        next_status_print = 100_000
        for i, (an_iter, bn_iter, metadata) in enumerate(self._iter_domains_with_cache(FIRST_STEP_MAX_ITERS)):
            try:
                key, _ = gcf_calculation_to_precision(an_iter, bn_iter, g_N_initial_key_length, FIRST_STEP_MIN_ITERS,
                                                      FIRST_STEP_BURST_NUMBER)
            except (ZeroInAn, NotConverging, ZeroDivisionError):
                continue

            if key in self.hash_table:  # find hits in hash table
                results.append(Match(key, metadata.an_coef, metadata.bn_coef))

            if i == next_status_print:  # print status.
                next_status_print += 1000
                print(
                    f'passed {i} out of {self.poly_domains.num_iterations} ' +
                    f'({round(100. * i / self.poly_domains.num_iterations, 2)}%). ' +
                    f' found so far {len(results)} results')
                print(f'currently at an = {metadata.an_coef} bn = {metadata.bn_coef}')

        if verbose:
            print(f'created results after {time() - start}s')
            print(f'found {len(results)} results')
        return results

    def _improve_results_precision(self, intermediate_results: List[Match], verbose=True):
        """
        For each results, calculate the GCD to a higher depth, and return the calculated result 
        with the original result.
        """
        precise_results = []
        counter = 0
        n_iterations = len(intermediate_results)
        key_factor = 10 ** g_N_verify_compare_length
        an_series_iter, bn_series_iter = self.poly_domains.get_calculation_method()

        for res in intermediate_results:
            counter += 1
            
            if counter % 10_000 == 0 and verbose:
                print('Calculated {} matches out of {} to a more precise value.'.format(
                    counter, n_iterations))

            an_iter = an_series_iter(res.rhs_an_poly, SECOND_STEP_MAX_ITERS, start_n=0)
            bn_iter = bn_series_iter(res.rhs_bn_poly, SECOND_STEP_MAX_ITERS, start_n=0)

            try:
                long_key, precision = gcf_calculation_to_precision(an_iter, bn_iter, g_N_verify_compare_length,
                                                                   min_iters=SECOND_STEP_MIN_ITERS,
                                                                   burst_number=SECOND_STEP_BURST_NUMBER)
            except NotConverging as e:
                print(f"{res} does not converge. Continuing...")
                print(e)
                continue
            except (ZeroInAn, ZeroDivisionError) as e:
                print(f" exception for {res}. Continuing...")
                print(e)
                continue

            rhs_val = mpmath.mpf(long_key) / key_factor
            rhs_str = mpmath.nstr(rhs_val, precision)

            precise_results.append((res, rhs_str, precision))

        print(f'finished calculating {n_iterations} to a more precise value.')
        return precise_results

    def _refine_results(self, precise_intermediate_results, verbose=True):
        """
        This step is identical to the one in Efficient enumerator. The two we're quite different in the
        rest of the implementation, so it didn't feel right to extend EfficientGCFEnumerator in this class.
        Hence, there is some code duplication here

        validate intermediate results to 100 digit precision
        :param precise_intermediate_results:  list of results from first enumeration
        :param verbose: if true print status.
        :return: final results.
        """
        results = []
        counter = 0
        n_iterations = len(precise_intermediate_results)

        for res, rhs_str, precision in precise_intermediate_results:
            counter += 1
            if (counter % 10_000 == 0 or counter == n_iterations) and verbose:
                print('Passed {} permutations out of {}. Found so far {} matches'.format(
                    counter, n_iterations, len(results)))

            try:
                all_matches = self.hash_table.evaluate(res.lhs_key)
            except KeyError:
                continue

            for i, match in enumerate(all_matches):
                # Trunc the LHS to the number of digits actually calculated.
                val_str = mpmath.nstr(match[0], precision + 1)[:-1]
                if val_str == rhs_str:
                    results.append(RefinedMatch(*res, i, match[1], match[2], precision))

        return results
