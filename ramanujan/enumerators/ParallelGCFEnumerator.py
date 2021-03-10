import itertools
import numpy as np
from time import time
from typing import List, Iterator, Callable

from ramanujan.constants import g_N_initial_search_terms
from .AbstractGCFEnumerator import Match, RefinedMatch
from .EfficientGCFEnumerator import EfficientGCFEnumerator


class ParallelGCFEnumerator(EfficientGCFEnumerator):
    """
    Parallel implementation of EfficientGCFEnumerator's _first_enumeration.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod  # Copied from EfficientGCFEnumerator, as it is a private method
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


    # Override
    def _first_enumeration(self, verbose: bool) -> List[Match]:
        """ See EfficientGCFEnumerator for documentation! """
        def efficient_gcf_calculation(shape: List[int], length: int) -> np.ndarray:
            """
            enclosure. a_, b_, and key_factor are used from outer scope.
            :param shape: common shape for arrays
            :param length: common dimension of a and b
            :return: keys for LHS hash table
            """
            prev_q = np.zeros(shape, dtype=np.float64)
            q = np.ones(shape, dtype=np.float64)
            next_q = np.empty(shape, dtype=np.float64)
            prev_p = np.ones(shape, dtype=np.float64)
            p = np.full(shape, a_[0, ..., np.newaxis], dtype=np.float64)  # newaxis for broadcasting
            next_p = np.empty(shape, dtype=np.float64)
            temp1 = np.empty(shape, dtype=np.float64)
            temp2 = np.empty(shape, dtype=np.float64)
            
            for i in range(1, length):
                temp1 = np.multiply(a_[i, ..., np.newaxis], q, out=temp1)
                temp2 = np.multiply(b_[np.newaxis, i], prev_q, out=temp2)
                next_q = np.add(temp1,  temp2, out=next_q)
                
                temp1 = np.multiply(a_[i, ..., np.newaxis], p, out=temp1)
                temp2 = np.multiply(b_[np.newaxis, i], prev_p, out=temp2)
                next_p = np.add(temp1, temp2, out=next_p)
                
                # Move references one step forward
                prev_q, q, next_q = q, next_q, prev_q
                prev_p, p, next_p = p, next_p, prev_p

            result = np.multiply(key_factor, p, out=temp1)
            result = np.divide(result, q, out=result, where=(q != 0))
            return np.trunc(result, out=result)
            
        start = time()
        key_factor = round(1 / self.threshold)
        counter = print_counter = 0  # number of permutations passed
        results = []  # list of intermediate results

        a_coef_list, an_list = self.__create_series_list(
            self.get_an_iterator(), self.create_an_series, filter_from_1=True)
        b_coef_list, bn_list = self.__create_series_list(
            self.get_bn_iterator(), self.create_bn_series, filter_from_1=True)
        
        num_iterations = len(an_list) * len(bn_list)
        if num_iterations == 0:
            print("Nothing to iterate over!")
            return []
    
        a_ = np.array(an_list, dtype=np.float64).T
        b_ = np.array(bn_list, dtype=np.float64).T
        
        if verbose:
            print(f'Created final enumerations filters after {time() - start:.2f}s')
            start = time()
            
        # calculate hash key of gcf value    
        many_keys = efficient_gcf_calculation((*a_.shape[1:], *b_.shape[1:]), a_.shape[0])
        print(type(many_keys[0,0]))
            
        if verbose:
            print(f"Calculations in {time() - start:.2f}s")
            start = time()
            
        for bind, b_coef in enumerate(b_coef_list):
            for aind, a_coef in enumerate(a_coef_list):
                key = int(many_keys[aind, bind])
                if key in self.hash_table:  # find hits in hash table (bottleneck)
                    results.append(Match(key, a_coef, b_coef))
            
            if verbose:
                counter += len(an_list)
                print_counter += len(an_list)
                if print_counter >= 1_000_000:  # print status.
                    print_counter = 0
                    prediction = (time() - start)*(num_iterations / counter)
                    time_left = (time() - start)*(num_iterations / counter - 1)
                    print(f"Passed {counter} out of {num_iterations} " +
                          f"({round(100. * counter / num_iterations, 2)}%). "
                          f"Found so far {len(results)} results. "
                          f"Time left ~{time_left:.0f}s of a total of {prediction:.0f}s")

        if verbose:
            print(f'created results after {time() - start:.2f}s')
        return results
            
