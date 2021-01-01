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
from utils import find_polynomial_series_coefficients
from LHSHashTable import LHSHashTable

# GCF enumerator type 
from RelativeGCFEnumerator import RelativeGCFEnumerator
from EfficentGCFEnumerator import EfficentGCFEnumerator
from MonomExponentialGCFEnumerator import MonomExponentialGCFEnumerator

enumerators = { 
    'relative': RelativeGCFEnumerator,
    'efficent': EfficentGCFEnumerator,
    'monoms': MonomExponentialGCFEnumerator
}

def multi_core_enumeration(sym_constant, lhs_search_limit, saved_hash, poly_a, poly_b, num_cores, splits_size,
                           create_an_series=None, create_bn_series=None, index=0, enumerator_type='efficent'):
    """
    function to run for each process. this also divides the work to tiles/
    :param sym_constant: sympy constant for search
    :param lhs_search_limit:  limit for hash table
    :param saved_hash: path to saved hash table
    :param poly_a: explained in docstring of __first_enumeration
    :param poly_b: explained in docstring of __first_enumeration
    :param num_cores: total number of cores used.
    :param splits_size: tile size for each process.
    we can think of the search domain as a n-d array with dim(poly_a) + dim(poly_b) dimensions.
    to split this efficiently we need the tile size. for each value in splits_size we take it as the tile size for a
    dimension of the search domain. for example, is split size is [4,5] then we will split the work in the first
    2 dimensions of the search domain to tiles of size [4,5].
    NOTICE - we do not verify that the tile size make sense to the number of cores used.
    :param index: index of core used.
    :param create_an_series: a custom function for creating a_n series with poly_a coefficients
    (default is create_series_from_compact_poly)
    :param create_bn_series: a custom function for creating b_n series with poly_b coefficients
    (default is create_series_from_compact_poly)
    :return: results
    """
    for s in range(len(splits_size)):
        if splits_size[s] == -1:  # no split
            continue
        if index == (num_cores - 1):  # last processor does more.
            poly_a[s] = poly_a[s][index * splits_size[s]:]
        else:
            try:
                poly_a[s] = poly_a[s][index * splits_size[s]:(index + 1) * splits_size[s]]
            except Exception as e:
                import ipdb
                ipdb.set_trace()

    import ipdb
    ipdb.set_trace()
    enumerator = enumerators[enumerator_type](sym_constant, lhs_search_limit, saved_hash, create_an_series, create_bn_series)
    results = enumerator.find_initial_hits(poly_a, poly_b, index == (num_cores - 1))
    return results


def multi_core_enumeration_wrapper(sym_constant, lhs_search_limit, poly_a, poly_b, num_cores, manual_splits_size=None,
                                   saved_hash=None, create_an_series=None, create_bn_series=None, enumerator_type='efficent'):
    """
    a wrapper for enumerating using multi processing.
    :param sym_constant: sympy constant for search
    :param lhs_search_limit: limit for hash table
    :param poly_a: explained in docstring of __first_enumeration
    :param poly_b: explained in docstring of __first_enumeration
    :param num_cores: total number of cores to be used.
    :param manual_splits_size: amount of work for each processor.
    manuals tiling (explained in docstring of multi_core_enumeration)
    by default we will split the work only along the first dimension. so the tile size will be
    [dim0 / n_cores, . , . , . , rest of dimensions].
    passing this manually can be useful for a large number of cores.
    :param saved_hash: path to saved hash table file if exists.
    :param create_an_series: a custom function for creating a_n series with poly_a coefficients
    (default is create_series_from_compact_poly)
    :param create_bn_series: a custom function for creating b_n series with poly_b coefficients
    (default is create_series_from_compact_poly)
    :return: results.
    """
    print(locals())
    if not os.path.isfile(saved_hash):
        if saved_hash is None:  # if no hash table given, build it here.
            saved_hash = 'tmp_hash.p'
        enumerator = enumerators[enumerator_type](sym_constant, lhs_search_limit, saved_hash)
        enumerator.hash_table.save()  # and save it to file (and global instance)
    else:
        if os.name != 'nt':  # if creation of process uses 'Copy On Write' we can benefit from it by
            # loading the hash table to memory here.
            enumerators[enumerator_type](sym_constant, lhs_search_limit, saved_hash)

    if manual_splits_size is None:  # naive work split
        manual_splits_size = [len(poly_a[0]) // num_cores]

    # built function for processes
    func = partial(multi_core_enumeration, sym_constant, lhs_search_limit, saved_hash, poly_a, poly_b, num_cores,
                   manual_splits_size, create_an_series, create_bn_series, enumerator_type=enumerator_type)

    if num_cores == 1:  # don't open child processes
        results = func(0)
    else:
        print(
            'starting Multi-Processor search.\n\tNOTICE- intermediate status prints will be done by processor 0 only.')
        pool = multiprocessing.Pool(num_cores)
        partial_results = pool.map(func, range(num_cores))
        results = []
        for r in partial_results:
            results += r

    enumerator = enumerators[enumerator_type](sym_constant, lhs_search_limit, saved_hash, create_an_series, create_bn_series)
    results = enumerator.refine_results(results)
    print(f'found {len(results)} results!')

    print('results in unicode:')
    enumerator.print_results(results, latex=False, convergence_rate=False)
    print('results in latex:')
    enumerator.print_results(results, latex=True, convergence_rate=False)

    results_in_latex = enumerator.convert_results_to_latex(results)
    generate_latex(file_name=f'results/{datetime.now().strftime("%m-%d-%Y--%H-%M-%S")}', eqns=results_in_latex)

    return results_in_latex
