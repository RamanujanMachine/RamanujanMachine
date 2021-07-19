import itertools
import mpmath
from typing import List, Iterator, Callable
from time import time

from ramanujan.utils.mobius import EfficientGCF
from ramanujan.constants import g_N_initial_search_terms, g_N_verify_terms, g_N_verify_compare_length
from .AbstractGCFEnumerator import AbstractGCFEnumerator, Match, RefinedMatch
from ramanujan.utils.utils import trunc_division


class EfficientGCFEnumerator(AbstractGCFEnumerator):
    """
    This enumerator maximizes on efficiency for calculating GCF to a precise dept.

    first enumeration will calculate the GCF to dept g_N_initial_search_terms (about 30),
    and compare against the lhs table to g_N_initial_key_length (usually 10 digits)
    results refining will calculate the GCF to dept g_N_verify_terms (usually 1000)
    and compare it with g_N_verify_compare_length (100) digits of the given expression
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def __create_series_list(coefficient_iter: Iterator,
                             series_generator: Callable[[List[int], int], List[int]],
                             filter_from_1=False) -> [List[int], List[int]]:
        coef_list = list(coefficient_iter)
        # create a_n and b_n series for coefficients.
        series_list = [series_generator(coef_list[i], g_N_initial_search_terms)
                       for i in range(len(coef_list))]
        # filter out all options resulting in '0' in any series term.
        if filter_from_1:
            series_filter = [0 not in an[1:] for an in series_list]
        else:
            series_filter = [0 not in an for an in series_list]
        series_list = list(itertools.compress(series_list, series_filter))
        coef_list = list(itertools.compress(coef_list, series_filter))
        return coef_list, series_list

    def _first_enumeration(self, verbose: bool):
        """
        This is usually the bottleneck of the search.
        We calculate general continued fractions of type K(bn,an).
        'an' and 'bn' are polynomial series.
        These polynomials take the form of n(n(..(n*c_1 + c_0) + c_2)..)+c_k.
        The polynomial families are supplied from self.poly_domains_generator

        For each an and bn pair, a gcf is calculated using efficient_gcf_calculation
        defined under this scope, and compared self.hash_tables for hits.

        :param verbose: if True print the status of calculation.
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
                q = a_[i] * q + b_[i] * prev_q
                p = a_[i] * p + b_[i] * prev_p
                prev_q = tmp_a
                prev_p = tmp_b
            # calculate hash key of gcf value
            return trunc_division(key_factor * p, q) if q != 0 else 0

        start = time()
        a_coef_iter = self.get_an_iterator()  # all coefficients possibilities for 'a_n'
        b_coef_iter = self.get_bn_iterator()
        size_b = self.get_bn_length()
        size_a = self.get_an_length()
        num_iterations = size_b * size_a
        key_factor = round(1 / self.threshold)

        counter = 0  # number of permutations passed
        print_counter = counter
        results = []  # list of intermediate results

        if size_a > size_b:  # cache {bn} in RAM, iterate over an
            b_coef_list, bn_list = self.__create_series_list(b_coef_iter, self.create_bn_series,
                                                             filter_from_1=True)
            real_bn_size = len(bn_list)
            num_iterations = (num_iterations // self.get_bn_length()) * real_bn_size
            if verbose:
                print(f'created final enumerations filters after {time() - start:.2f}s')
            start = time()
            for a_coef in a_coef_iter:
                an = self.create_an_series(a_coef, g_N_initial_search_terms)
                if 0 in an[1:]:  # a_0 is allowed to be 0.
                    counter += real_bn_size
                    print_counter += real_bn_size
                    continue
                for bn_coef in zip(bn_list, b_coef_list):
                    a_ = an
                    b_ = bn_coef[0]
                    key = efficient_gcf_calculation()  # calculate hash key of gcf value

                    if key in self.hash_table:  # find hits in hash table
                        results.append(Match(key, a_coef, bn_coef[1]))
                    if verbose:
                        counter += 1
                        print_counter += 1
                        if print_counter >= 1_000_000:  # print status.
                            print_counter = 0
                            prediction = (time() - start)*(num_iterations / counter)
                            time_left = (time() - start)*(num_iterations / counter - 1)
                            print(f'Passed {counter} out of {num_iterations} '
                                  f'({round(100. * counter / num_iterations, 2)}%). '
                                  f'Found so far {len(results)} results. \n'
                                  f'Time left ~{time_left:.0f}s of a total of {prediction:.0f}s')

        else:  # cache {an} in RAM, iterate over bn
            a_coef_list, an_list = self.__create_series_list(a_coef_iter, self.create_an_series,
                                                             filter_from_1=True)
            real_an_size = len(an_list)
            num_iterations = (num_iterations // self.get_an_length()) * real_an_size
            if verbose:
                print(f'created final enumerations filters after {time() - start:.2f}s')
            start = time()
            for b_coef in b_coef_iter:
                bn = self.create_bn_series(b_coef, g_N_initial_search_terms)
                if 0 in bn[1:]:
                    counter += real_an_size
                    print_counter += real_an_size
                    continue
                for an_coef in zip(an_list, a_coef_list):
                    a_ = an_coef[0]
                    b_ = bn
                    key = efficient_gcf_calculation()  # calculate hash key of gcf value

                    if key in self.hash_table:  # find hits in hash table
                        results.append(Match(key, an_coef[1], b_coef))
                    if verbose:
                        counter += 1
                        print_counter += 1
                        if print_counter >= 1_000_000:  # print status.
                            print_counter = 0
                            prediction = (time() - start)*(num_iterations / counter)
                            time_left = (time() - start)*(num_iterations / counter - 1)
                            print(f'Passed {counter} out of {num_iterations} '
                                  f'({round(100. * counter / num_iterations, 2)}%). '
                                  f'Found so far {len(results)} results. \n'
                                  f'Time left ~{time_left:.0f}s of a total of {prediction:.0f}s')

        if verbose:
            print(f'created results after {time() - start:.2f}s')
        return results

    def _improve_results_precision(self, intermediate_results: List[Match], verbose=True):
        """
        For each results, calculate the GCF to a higher dept, and return the calculated result 
        with the original result.
        """
        precise_results = []
        counter = 0
        n_iterations = len(intermediate_results)
        for res in intermediate_results:
            counter += 1
            if (counter % 10_000 == 0 or counter == n_iterations) and verbose:
                print('Calculated {} matches out of {} to a more precise value.'.format(
                    counter, n_iterations))

            # create a_n, b_n with huge length, calculate gcf, and verify result.
            an = self.create_an_series(res.rhs_an_poly, g_N_verify_terms)
            bn = self.create_bn_series(res.rhs_bn_poly, g_N_verify_terms)
            gcf = EfficientGCF(an, bn)
            rhs_str = mpmath.nstr(gcf.evaluate(), g_N_verify_compare_length)

            precise_results.append((res, rhs_str))

        return precise_results

    def _refine_results(self, precise_intermediate_results, verbose=True):
        """
        validate intermediate results to 100 digit precision
        :param precise_intermediate_results:  list of results from first enumeration
        :param verbose: if true print status.
        :return: final results.
        """
        results = []
        counter = 0
        n_iterations = len(precise_intermediate_results)
        constant_vals = [const() for const in self.constants_generator]

        for res, rhs_str in precise_intermediate_results:
            counter += 1
            if (counter % 10_000 == 0 or counter == n_iterations) and verbose:
                print('Passed {} permutations out of {}. Found so far {} matches'.format(
                    counter, n_iterations, len(results)))
            try:
                all_matches = self.hash_table.evaluate(res.lhs_key)
                # check if all values encountered are not inf or nan
                if not all([not (mpmath.isinf(val) or mpmath.isnan(val))
                            for val, _, _ in all_matches]):  # safety
                    print('Something wicked happened!')
                    print(f'Encountered a NAN or inf in LHS db, at {res.lhs_key}, {constant_vals}')
                    continue
            except (ZeroDivisionError, KeyError):
                # if there was an exception here, there is no need to halt the entire execution,
                # but only note it to the user
                continue

            for i, match in enumerate(all_matches):
                val_str = mpmath.nstr(match[0], g_N_verify_compare_length)
                if val_str == rhs_str:
                    # This patch is meant to allow support for multiple matches for an
                    # LHS key, it will later be used to determine which item in the LHS dict
                    # was matched
                    results.append(RefinedMatch(*res, i, match[1], match[2]))

        return results
