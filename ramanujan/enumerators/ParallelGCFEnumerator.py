import itertools
import mpmath as mp
import numpy as np
from time import time
from typing import List, Iterator, Callable, Optional

from ramanujan.constants import g_N_initial_search_terms
from .AbstractGCFEnumerator import Match, RefinedMatch
from .EfficientGCFEnumerator import EfficientGCFEnumerator


MAX_RAM = 0.1 * 2.**30  # 0.1 GB


def calculate_RAM_usage(shape):
    """ 
    Calculate a lower bound to the RAM needed to allocate the arrays for efficient_gcf_calculation
    
    8 arrays of shape full of 64-bit floats. 
    There are also some lists whose space is hard to estimate. 
    
    :param shape: see efficient_gcf_calculation
    :param length: see efficient_gcf_calculation
    :return float: RAM usage in bytes
    """
    c = 200.
    return 8. * mp.fprod(shape) * 8. + c * mp.fsum(shape) * 8.


class ParallelGCFEnumerator(EfficientGCFEnumerator):
    """
    Parallel implementation of EfficientGCFEnumerator's _first_enumeration.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def __create_series_list(coefficient_iter: Iterator,
                             series_generator: Callable[[List[int], int], List[int]],
                             filter_from_1=False,
                             iterations: Optional[int]=None) -> [List[int], List[int]]:
        """ Generate coefficients and respective series """
        coef_list = list()
        series_list = list()
        # create a_n and b_n series for coefficients.
        for i, coef in enumerate(coefficient_iter):
            if i == iterations:
                break
            an = series_generator(coef, g_N_initial_search_terms)
            # filter out all options resulting in '0' in any series term.
            series_filter = 0 not in an[1:] if filter_from_1 else 0 not in an
            if series_filter:
                coef_list.append(coef)
                series_list.append(an)
            
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
            p = np.full(shape, a_[0, ..., np.newaxis], dtype=np.float64)
            next_p = np.empty(shape, dtype=np.float64)
            temp1 = np.empty(shape, dtype=np.float64)
            temp2 = np.empty(shape, dtype=np.float64)
            # new axis for numpy broadcasting
            
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
        counter = 0  # number of permutations passed
        print_counter = calc_time = chunks_done = 0
        results = []  # list of intermediate results

        asize = self.get_an_length()
        bsize = self.get_bn_length()
        
        num_iterations = asize * bsize
        if num_iterations == 0:
            print("Nothing to iterate over!")
            return []
    
        # Split task into chunks
        min_chunks = round(np.ceil(calculate_RAM_usage((asize, bsize)) / MAX_RAM))
        if min_chunks < max(asize, bsize):  # Iterate over intervals on the longer axis
            achunk = asize if asize < bsize else np.int(np.ceil(asize / min_chunks))
            bchunk = bsize if asize >= bsize else np.int(np.ceil(bsize / min_chunks))
        else:  # Iterate over intervals on the longer axis for each on the shorter axis
            achunk = 1 if asize < bsize else np.int(np.ceil(asize * bsize / min_chunks))
            bchunk = 1 if asize >= bsize else np.int(np.ceil(asize * bsize / min_chunks))
        
        if verbose:
            chunks_total = round(np.ceil(asize / achunk) * np.ceil(bsize / bchunk))
            print(f'Created final enumerations filters after {time() - start:.2f}s')
            print(f"Doing {num_iterations} searches in up to {chunks_total} chunks. "
                  f"This might take some time. ")
            start_results = time()    
        
        if asize < bsize:
            small_chunk, large_chunk = achunk, bchunk
            small_size, large_size = asize, bsize
            small_iterator, large_iterator = self.get_an_iterator, self.get_bn_iterator
            small_series, large_series = self.create_an_series, self.create_bn_series
            small_poly, large_poly = {"coef": 0, "series": 0}, {"coef": 0, "series": 0}
            a_poly, b_poly = small_poly, large_poly  # Link the dictionaries
        else:
            small_chunk, large_chunk = bchunk, achunk
            small_size, large_size = bsize, asize
            small_iterator, large_iterator = self.get_bn_iterator, self.get_an_iterator
            small_series, large_series = self.create_bn_series, self.create_an_series
            small_poly, large_poly = {"coef": 0, "series": 0}, {"coef": 0, "series": 0}
            b_poly, a_poly = small_poly, large_poly  # Link the dictionaries
        
        # Compute matches
        large_iter = large_iterator()
        for _ in range(0, large_size, large_chunk):
            large_poly["coef"], large_poly["series"] = self.__create_series_list(
                    large_iter, large_series, filter_from_1=True, iterations=large_chunk)
            if len(large_poly["series"]) == 0:  # exhausted or all include 0
                continue
            
            small_iter = small_iterator()
            for _ in range(0, small_size, small_chunk):
                start = time()
                
                small_poly["coef"], small_poly["series"] = self.__create_series_list(
                    small_iter, small_series, filter_from_1=True, iterations=small_chunk)
                if len(small_poly["series"]) == 0:  # exhausted or all include 0
                    continue 
                
                a_ = np.array(a_poly["series"], dtype=np.float64).T
                b_ = np.array(b_poly["series"], dtype=np.float64).T
                shape = (a_.shape[1], b_.shape[1])
                
                # calculate hash key of gcf value  
                many_keys = efficient_gcf_calculation(shape, a_.shape[0])
                    
                if verbose:
                    calc_time += time() - start
                    print(f"Calculations in {time() - start:.2f}s")
                    start = time()
                    print_counter = 0  # Less prints with many chunks
                    chunks_done += 1
                    
                for aind in range(shape[0]):
                    for bind in range(shape[1]):
                        key = int(many_keys[aind, bind])
                        if key in self.hash_table:  # find hits in hash table (bottleneck)
                            results.append(Match(key, a_poly["coef"][aind], b_poly["coef"][bind]))
                    
                    if verbose:
                        counter += shape[1]
                        print_counter += shape[1]
                        if print_counter >= 1_000_000:  # print status.
                            print_counter = 0
                            prediction = ((time() - start_results - calc_time)
                                            *(num_iterations / counter)
                                            + calc_time * chunks_total / chunks_done)
                            time_left = ((time() - start_results - calc_time)
                                            *(num_iterations / counter - 1)
                                            + calc_time * (chunks_total / chunks_done - 1))
                            print(f"Passed {counter:n} out of {num_iterations:n} "
                                  f"({round(100. * counter / num_iterations, 2)}%). "
                                  f"Found so far {len(results)} results. \n"
                                  f"Time left ~{time_left:.0f}s of a total of {prediction:.0f}s")
                
                if verbose: # Chunk complete
                    prediction = (time() - start_results)*(num_iterations / counter)
                    if prediction < 120:
                        prediction = f"{prediction:.0f}s"
                    elif prediction < 60*60:
                        prediction = f"{prediction//60:.0f}min {prediction%60:.0f}s"
                    else:
                        prediction = f"{prediction//3600:.0f}h {prediction%3600//60:.0f}min"
                    
                    time_left = (time() - start_results)*(num_iterations / counter - 1)
                    if time_left < 120:
                        time_left = f"{time_left:.0f}s"
                    elif time_left < 60*60:
                        time_left = f"{time_left//60:.0f}min {time_left%60:.0f}s"
                    else:
                        time_left = f"{time_left//3600:.0f}h {time_left%3600//60:.0f}min"
                    
                    print(f"Chunk {chunks_done}/{chunks_total} done. "
                          f"({round(100. * chunks_done / chunks_total, 2)}%). "
                          f"Time left {time_left} of a total of {prediction}")
                            

        if verbose:
            print(f'created results after {time() - start_results:.2f}s')
        return results
            
