import pickle
import mpmath
import numpy as np
from sympy import lambdify
import sympy
from time import time
import itertools
from massey import create_series_from_compact_poly
from mobius import GeneralizedContinuedFraction, MobiusTransform, EfficientGCF
import os
from collections import namedtuple
from typing import List
from math import gcd
from typing import TypeVar, Iterator
import multiprocessing
from functools import partial

# intermediate result - coefficients of lhs transformation, and compact polynomials for seeding an and bn series.
Match = namedtuple('Match', 'lhs_coefs rhs_an_poly rhs_bn_poly')


class GlobalHashTableInstance:
    def __init__(self):
        """
        python processes don't share memory. so when using multiprocessing the hash table will be duplicated.
        to try and avoid this, we initiate a global instance of the hash table.
        hopefully this will useful when running on linux (taking advantage of Copy On Write).
        this has not been tested yet on linux.
        on windows it has no effect when multiprocessing.
        """
        self.hash = {}
        self.name = ''


# global instance
hash_instance = GlobalHashTableInstance()


class LHSHashTable(object):

    def __init__(self, search_range_top, search_range_bottom, const_val, threshold) -> None:
        """
        hash table for LHS. storing values in the form of (ax + b)/(cx + d)
        :param search_range_top: range for values a,b.
        :param search_range_bottom: range for value c,d.
        :param const_val: constant for x.
        :param threshold: decimal threshold for comparison. in fact, the keys for hashing will be the first
                            -log_{10}(threshold) digits of the value. for example, if threshold is 1e-10 - then the
                            first 10 digits will be used as the hash key.
        """
        self.s = {}
        self.threshold = threshold
        for a in search_range_top:
            for b in search_range_top:
                for c in search_range_bottom:
                    for d in search_range_bottom:
                        if gcd(gcd(a, b), gcd(c, d)) != 1:  # don't store values that already exist
                            continue
                        denominator = c*const_val + d
                        numerator = a*const_val + b
                        if denominator == 0 or numerator == 0:  # don't store nan or 0.
                            continue
                        val = numerator / denominator
                        if mpmath.isnan(val) or mpmath.isinf(val):  # safety check
                            continue
                        if ((c + d) != 0) and mpmath.almosteq(val, ((mpmath.mpf(a)+mpmath.mpf(b))/(c+d))):
                            # don't store values that are independent of the constant (e.g. rational numbers)
                            continue
                        key = int(val / self.threshold)
                        if key in self.s:
                            continue
                        self.s[key] = np.array([[a, b], [c, d]], dtype=object)  # store key and transformation

    def __contains__(self, item):
        """
        operator 'in'
        :param item: key
        :return: true of false
        """
        return item in self.s

    def __getitem__(self, item):
        """
        operator []
        :param item: key
        :return: transformation of x
        """
        return self.s[item]

    def __eq__(self, other):
        """
        operator ==
        :param other: other hash table.
        :return:
        """
        if type(other) != type(self):
            return False
        ret = self.threshold == other.threshold
        ret &= sorted(self.s.keys()) == sorted(other.s.keys())
        return ret

    def save(self, name):
        """
        save the hash table as file
        :param name: path for file.
        """
        if hash_instance.name != name:  # save to global instance.
            hash_instance.hash = self
            hash_instance.name = name
        with open(name, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load_from(cls, name):
        """
        load hash table from file (or global instance)
        :param name:
        :return:
        """
        if hash_instance.name == name:
            print('loading instance')
            return hash_instance.hash   # hopefully on linux this will not make a copy.
        else:
            with open(name, 'rb') as f:
                print('not loading instance')
                ret = pickle.load(f)
                hash_instance.hash = ret    # save in instance
                hash_instance.name = name
        return ret


T = TypeVar('T')    # template


def chunks(iterator: Iterator[T], n: int) -> Iterator[Iterator[T]]:
    """
    creating chunk iterator
    ref - https://dev.to/orenovadia/solution-chunked-iterator-python-riddle-3ple
    """
    for first in iterator:  # take one item out (exits loop if `iterator` is empty)
        rest_of_chunk = itertools.islice(iterator, 0, n - 1)
        yield itertools.chain([first], rest_of_chunk)  # concatenate the first item back


class EnumerateOverGCF(object):
    def __init__(self, sym_constant, lhs_search_limit, saved_hash=''):
        """
        initialize search engine.
        basically, this is a 3 step procedure:
        1) load / initialize lhs hash table.
        2) first enumeration - enumerate over all rhs combinations, find hits in lhs hash table.
        3) refine results - take results from (2) and validate them to 100 decimal digits.
        :param sym_constant: sympy constant
        :param lhs_search_limit: range of coefficients for left hand side.
        :param saved_hash: path to saved hash.
        """
        self.threshold = 1e-10  # key length
        self.enum_dps = 50      # working decimal precision for first enumeration
        self.verify_dps = 2000   # working decimal precision for validating results.
        self.lhs_limit = lhs_search_limit
        self.const_sym = sym_constant
        self.const_val = lambdify((), sym_constant, modules="mpmath")
        if saved_hash == '':
            print('no previous hash table given, initializing hash table...')
            with mpmath.workdps(self.enum_dps):
                start = time()
                self.hash_table = LHSHashTable(
                    range(self.lhs_limit+1),    # a,b range (allow only non-negative)
                    range(-self.lhs_limit, self.lhs_limit+1),   # c,d range
                    self.const_val(),   # constant
                    self.threshold)     # length of key
                end = time()
                print('that took {}s'.format(end-start))
        else:
            self.hash_table = LHSHashTable.load_from(saved_hash)

    def __first_enumeration(self, poly_a: List[List], poly_b: List[List], print_results: bool):
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
        start = time()
        neg_poly_b = [[-i for i in b] for b in poly_b]      # for b_n include negative terms
        a_coef_list = list(itertools.product(*poly_a))      # all coefficients possibilities for 'a_n'
        b_coef_list = list(itertools.product(*poly_b)) + list(itertools.product(*neg_poly_b))
        num_iterations = len(a_coef_list) * len(b_coef_list)  # total number of permutations

        # create a_n and b_n series fro coefficients.
        an_list = [create_series_from_compact_poly(a_coef_list[i], 32) for i in range(len(a_coef_list))]
        bn_list = [create_series_from_compact_poly(b_coef_list[i], 32) for i in range(len(b_coef_list))]
        if print_results:
            print('created series after {} s'.format(time() - start))

        # filter out all options resulting in '0' in any series term.
        an_filter = [0 not in an for an in an_list]
        bn_filter = [0 not in bn for bn in bn_list]
        an_list = list(itertools.compress(an_list, an_filter))
        a_coef_list = list(itertools.compress(a_coef_list, an_filter))
        bn_list = list(itertools.compress(bn_list, bn_filter))
        b_coef_list = list(itertools.compress(b_coef_list, bn_filter))

        # create another product of a_n options and b_n options
        an_bn_list = itertools.product(*[an_list, bn_list])
        a_b_coefs = itertools.product(*[a_coef_list, b_coef_list])
        if print_results:
            print('created final enumerations filters after {} s'.format(time() - start))

        counter = 0     # number of permutations passed
        results = []    # list of intermediate results
        batch_size = 10000  # chosen empirically, balance memory usage and cpu performance.
        start = time()
        for an_bn in zip(chunks(an_bn_list, batch_size), chunks(a_b_coefs, batch_size)):    # enumerate in batches.
            gcf_list = [EfficientGCF(an, bn) for (an, bn) in an_bn[0]]               # create gcf from a_n and b_n
            keys = [int(gcf.evaluate() / self.threshold) for gcf in gcf_list]        # calculate hash key of gcf value
            hits = [k in self.hash_table for k in keys]                              # find hits in hash table
            key_hits = list(itertools.compress(keys, hits))                          # filter by hits.
            a_b_coefs_hits = list(itertools.compress(an_bn[1], hits))                # and save them.
            results.extend([Match(self.hash_table[m[0]], m[1][0], m[1][1]) for m in zip(key_hits, a_b_coefs_hits)])
            if print_results:
                counter += batch_size
                if counter % 100000 == 0:   # print status.
                    print('passed {} out of {}. found so far {} results'.format(counter, num_iterations, len(results)))
        if print_results:
            print('created results after {} s'.format(time() - start))
        return results

    def __refine_results(self, intermediate_results: List[Match], print_results=True):
        """
        validate intermediate results to 100 digit precision
        :param intermediate_results:  list of results from first enumeration
        :param print_results: if true print status.
        :return: final results.
        """
        results = []
        counter = 0
        n_iterations = len(intermediate_results)
        for r in intermediate_results:
            counter += 1
            if (counter % 10) == 0 and print_results:
                print('passed {} permutations out of {}. found so far {} matches'.format(
                    counter, n_iterations, len(results)))
            t = MobiusTransform(r.lhs_coefs)
            try:
                val = t(self.const_val())
                if mpmath.isinf(val) or mpmath.isnan(val):  # safety
                    continue
                if mpmath.almosteq(val, t(1), 1/(self.verify_dps//20)):
                    # don't keep results that are independent of the constant
                    continue
            except ZeroDivisionError:
                continue

            # create a_n, b_n with huge length, calculate gcf, and verify result.
            an = create_series_from_compact_poly(r.rhs_an_poly, 1000)
            bn = create_series_from_compact_poly(r.rhs_bn_poly, 1000)
            gcf = EfficientGCF(an, bn)
            val_str = mpmath.nstr(val, 100)
            rhs_str = mpmath.nstr(gcf.evaluate(), 100)
            if val_str == rhs_str:
                results.append(r)
        return results

    def print_results(self, results: List[Match], latex=False):
        """
        pretty print the the results.
        :param results: list of final results as received from refine_results.
        :param latex: if True print in latex form, otherwise pretty print in unicode.
        """
        for r in results:
            an = create_series_from_compact_poly(r.rhs_an_poly, 1000)
            bn = create_series_from_compact_poly(r.rhs_bn_poly, 1000)
            print_length = max(max(len(r.rhs_an_poly), len(r.rhs_bn_poly)), 5)
            gcf = GeneralizedContinuedFraction(an, bn)
            t = MobiusTransform(r.lhs_coefs)
            sym_lhs = sympy.simplify(t.sym_expression(self.const_sym))
            if not latex:
                print('lhs: ')
                sympy.pprint(sym_lhs)
                print('rhs: ')
                gcf.print(print_length)
                print('lhs value: ' + mpmath.nstr(t(self.const_val()), 50))
                print('rhs value: ' + mpmath.nstr(gcf.evaluate(), 50))
            else:
                result = sympy.Eq(sym_lhs, gcf.sym_expression(print_length))
                print('$$ ' + sympy.latex(result) + ' $$')

    def find_hits(self, poly_a: List[List], poly_b: List[List], print_results=True):
        """
        use search engine to find results (steps (2) and (3) explained in __init__ docstring)
        :param poly_a: explained in docstring of __first_enumeration
        :param poly_b: explained in docstring of __first_enumeration
        :param print_results: if true, pretty print results at the end.
        :return: final results.
        """
        with mpmath.workdps(self.enum_dps):
            if print_results:
                print('starting preliminary search...')
            start = time()
            results = self.__first_enumeration(poly_a, poly_b, print_results)   # step (2)
            end = time()
            if print_results:
                print('that took {}s'.format(end-start))
        with mpmath.workdps(self.verify_dps):
            if print_results:
                print('starting to verify results...')
            start = time()
            refined_results = self.__refine_results(results, print_results)  # step (3)
            end = time()
            if print_results:
                print('that took {}s'.format(end-start))
            if print_results:
                self.print_results(refined_results)
        return refined_results


def multi_core_enumeration(sym_constant, lhs_search_limit, saved_hash, poly_a, poly_b, num_cores, splits_size, index):
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
    :return: results
    """
    for s in range(len(splits_size)):
        if index == (num_cores-1):
            poly_a[s] = poly_a[s][index * splits_size[s]:]
        else:
            poly_a[s] = poly_a[s][index * splits_size[s]:(index+1) * splits_size[s]]
    enumerator = EnumerateOverGCF(sym_constant, lhs_search_limit, saved_hash)
    results = enumerator.find_hits(poly_a, poly_b, index == 0)
    enumerator.print_results(results, True)
    return results


def multi_core_enumeration_wrapper(sym_constant, lhs_search_limit, poly_a, poly_b, num_cores, manual_splits_size=None,
                                   saved_hash=None):
    """

    :param sym_constant: sympy constant for search
    :param lhs_search_limit: limit for hash table
    :param poly_a: explained in docstring of __first_enumeration
    :param poly_b: explained in docstring of __first_enumeration
    :param num_cores: total number of cores to be used.
    :param manual_splits_size: manuals tiling (explained in docstring of multi_core_enumeration)
    by default we will split the work only along the first dimension. so the tile size will be
    [dim0 / n_cores, . , . , . , rest of dimensions].
    passing this manually can be useful for a large number of cores.
    :param saved_hash: path to saved hash table file if exists.
    :return: results.
    """
    if (saved_hash is None) or (not os.path.isfile(saved_hash)):
        if saved_hash is None:  # if no hash table given, build it here.
            saved_hash = 'tmp_hash.p'
        enumerator = EnumerateOverGCF(sym_constant, lhs_search_limit)
        enumerator.hash_table.save(saved_hash)   # and save it to file (and global instance)
    else:
        if os.name != 'nt':     # if creation of process uses 'Copy On Write' we can benefit from it by
            # loading the hash table to memory here.
            EnumerateOverGCF(sym_constant, lhs_search_limit, saved_hash)

    if manual_splits_size is None:  # naive work split
        manual_splits_size = [len(poly_a[0]) // num_cores]

    # built function for processes
    func = partial(multi_core_enumeration, sym_constant, lhs_search_limit, saved_hash, poly_a, poly_b, num_cores,
                   manual_splits_size)

    if num_cores == 1:  # don't open child processes
        results = func(0)
        print('found {} results!'.format(len(results)))
    else:
        pool = multiprocessing.Pool(num_cores)
        results = pool.map(func, range(num_cores))
        print('found {} results!'.format(sum([len(results[i]) for i in range(num_cores)])))
    return results


# TODO - create api for this.
if __name__ == "__main__":
    final_results = multi_core_enumeration_wrapper(sympy.zeta(2),     # constant to run on
                                                   20,          # lhs limit
                                                   [[i for i in range(16)]]*3,  # a_n polynomial coefficients
                                                   [[i for i in range(10)]]*5,  # b_n polynomial coefficients
                                                   4,           # number of cores to run on
                                                   None,        # use naive tiling
                                                   os.path.join('hash_tables', 'zeta2_20_hash.p')  # existing hash table
                                                   )

    #with open('results_of_e_30', 'wb') as file:
    #    pickle.dump(final_results, file)
