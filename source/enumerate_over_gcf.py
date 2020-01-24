import pickle
import mpmath
import numpy as np
from sympy import lambdify, floor
import sympy
from time import time
import itertools
from series_generators import create_series_from_compact_poly, create_zeta_bn_series, create_series_from_shift_reg
from massey import slow_massey
from mobius import GeneralizedContinuedFraction, MobiusTransform, EfficientGCF
import os
from collections import namedtuple
from typing import List
from math import gcd
from typing import TypeVar, Iterator
import multiprocessing
from functools import partial
from constants import redundant_cycles

# intermediate result - coefficients of lhs transformation, and compact polynomials for seeding an and bn series.
Match = namedtuple('Match', 'lhs_coefs rhs_an_poly rhs_bn_poly')

def clear_end_zeros(items):

    while items[-1] == 0:
        items.pop()


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
                        denominator = c * const_val + d
                        numerator = a * const_val + b
                        if denominator == 0 or numerator == 0:  # don't store nan or 0.
                            continue
                        val = numerator / denominator
                        if mpmath.isnan(val) or mpmath.isinf(val):  # safety check
                            continue
                        if ((c + d) != 0) and mpmath.almosteq(val, ((mpmath.mpf(a) + mpmath.mpf(b)) / (c + d))):
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
            return hash_instance.hash  # hopefully on linux this will not make a copy.
        else:
            with open(name, 'rb') as f:
                print('not loading instance')
                ret = pickle.load(f)
                hash_instance.hash = ret  # save in instance
                hash_instance.name = name
        return ret


T = TypeVar('T')  # template


def chunks(iterator: Iterator[T], n: int) -> Iterator[Iterator[T]]:
    """
    creating chunk iterator
    ref - https://dev.to/orenovadia/solution-chunked-iterator-python-riddle-3ple
    """
    for first in iterator:  # take one item out (exits loop if `iterator` is empty)
        rest_of_chunk = itertools.islice(iterator, 0, n - 1)
        yield itertools.chain([first], rest_of_chunk)  # concatenate the first item back


class SignedRcfEnumeration(object):
    """
            Initialize search engine.
            Basically, this is a 3 step procedure:
            1) Enumerates LHS symbolic expressions of rational functions of the constant, and non repeating sign periods.
            2) Iterates through domain. With low precision extracts a series. Checks if massey-pretty. Saves hits.
            3) Refine results - takes results from (2) and validate them to 100 decimal digits.
            Note that the structure of the enumeraion (and in fact of the problem), makes it possible to divide
            any domain to separate domains with respect to the signed cycles, but not to the
            :param sym_constant: sympy constant
            :param coefficients_limit: range of coefficients for the rational function on the LHS.
            :param cycle_len_range: range of lengths for the sign sequence's period.
            :param depth: Number of elements of a series to extract. Relates to length of typical LFSRs  of the consant.
            :param min_deg: Used to exclude lower degree polynomials.
            :param prime: Prime number in use by Massey algorithm.
            """
    def __init__(self, sym_constant, coefficients_limit, cycle_len_range,  depth, poly_deg, min_deg=None, prime=199):
        self.beauty_standard = 15
        self.enum_dps = 300
        self.verify_dps = 1000
        self.coeff_lim = coefficients_limit
        self.min_cycle_len = cycle_len_range[0]
        self.max_cycle_len = cycle_len_range[1]
        self.poly_deg = poly_deg
        self.min_deg = min_deg
        self.const_sym = sym_constant
        self.const_val = lambdify((), sym_constant, modules="mpmath")
        self.depth = depth
        self.verify_depth = 1000
        self.prime = prime


    def create_sign_seq_enumeration(self):
        """
        Creates a list of all possible sign sequences.
        Uses either a pre defined dictionary or a function to skip redundancies.
        """
        sign_seqs = []
        for cyc_len in range(self.min_cycle_len, self.max_cycle_len+1):
            sign_seqs = sign_seqs + list(itertools.product([-1,1],repeat=cyc_len))
        sign_seqs = [seq for seq in sign_seqs if not redundant_cycles[''.join([str(c) for c in seq])]]
        return sign_seqs

    def create_rational_symbol(self, numerator, denominator):
        numer_deg = len(numerator) -1
        denom_deg = len(denominator) -1
        numer_sym = 0
        denom_sym = 0
        for i in range(numer_deg + 1):
            numer_sym += numerator[i]*(self.const_sym**i)
        for i in range(denom_deg + 1):
            denom_sym += denominator[i]*(self.const_sym**i)

        return numer_sym/denom_sym


    def create_rational_variations_enum(self):
        """
        Creates a list of all possible rational expressions for the LHS.
        Expressions saved as sympy-simplified, positive expressions to reduce redundancy.
        Additional checks are performed to exclude degenerated cases.
        """
        print("Starting enumeration over LHS")
        start = time()
        coeffs = [i for i in range(-self.coeff_lim, self.coeff_lim+1)]
        if self.min_deg is not None:
            numerators = [list(numer) for numer in list(itertools.product(coeffs, repeat=self.poly_deg + 1)) \
                          if len(numer) >= self.min_deg+1]
            denominators = [list(denom) for denom in list(itertools.product(coeffs, repeat=self.poly_deg + 1)) \
                          if len(denom) >= self.min_deg+1]
        else:
            numerators = [list(numer) for numer in list(itertools.product(coeffs, repeat=self.poly_deg+1))]
            denominators = [list(denom) for denom in list(itertools.product(coeffs, repeat=self.poly_deg+1))]
        variations = itertools.product(numerators,denominators)
        expressions = set()
        for var in variations:
            numer = var[0]
            denom = var[1]
            if denom == [0 for i in denom] or numer == [0 for i in numer]:
                continue
            var_sym = self.create_rational_symbol(numer, denom)
            var_sym = sympy.simplify(var_sym)
            if var_sym == floor(var_sym): # degenerate integer.
                continue
            if abs(var_sym) not in expressions:
                expressions.add(abs(var_sym))
        print("Finished enumerations. Took {}  seconds".format(round(time()-start,2)))
        return expressions


    def find_signed_rcf_conj(self):
        """
        Builds the final domain.
        Iterates throgh the domain:
        extraction->massey->check->save.
        Additional checks are performed to exclude degenerated cases.
        """
        inter_results = []
        rational_variations = self.create_rational_variations_enum()
        sign_seqs = list(self.create_sign_seq_enumeration())
        domain_size = len(rational_variations) * len(sign_seqs)
        print("De-Facto Domain Size is: {}".format(domain_size))
        two_pc = max(domain_size // 20, 5)
        cnt = 0
        start = time()
        #iterate
        for instance in itertools.product(rational_variations, sign_seqs):
            cnt += 1
            var, sign_cyc = instance[0], list(instance[1])
            var_gen = lambdify((), var, modules="mpmath")
            seq_len = len(sign_cyc)
            if cnt % two_pc == 0:
                print("{}% of domain searched.\n".format(round(100 * cnt / domain_size, 2)))
                print("{} minutes passed.\n".format(round((time() - start) / 60, 2)))
            b_ = (sign_cyc * (self.depth//seq_len))[:self.depth]
            with mpmath.workdps(self.enum_dps):
                signed_rcf = GeneralizedContinuedFraction.from_irrational_constant(const_gen=var_gen, b_=b_)
            a_ = signed_rcf.a_
            if 0 in a_:
                continue
            if len(a_) < self.depth:
                continue
            a_sr = list(slow_massey(a_, self.prime))
            clear_end_zeros(a_sr)
            if len(a_sr) < self.beauty_standard:
                inter_results.append([var, sign_cyc, a_, a_sr])
        return inter_results


    def verify_results(self, results):
        """
        Validate intermediate results to 100 digit precision
        If a numeric value appears multiple times, the first is kept as valid. The rest saved as duplicates for later.
        """
        verified = []
        duplicates = {}
        res_set = set()
        for res in results:
            var_gen = lambdify((), res[0], modules="mpmath")
            a_ = create_series_from_shift_reg(res[3], res[2][:(len(res[3])-1)], self.verify_depth)
            b_ = (res[1] * (self.verify_depth // len(res[1])))[:self.verify_depth]
            gcf = EfficientGCF(a_, b_)
            with mpmath.workdps(self.verify_dps):
                lhs_str = mpmath.nstr(var_gen(), 100)
                rhs_val = gcf.evaluate()
                rhs_str = mpmath.nstr(rhs_val, 100)
                if rhs_str != lhs_str:
                    continue
                key = lhs_str
            if key not in res_set:
                res_set.add(key)
                verified.append(res)
                duplicates[key] = []
            else:
                duplicates[key].append(res)
                continue
        return verified, duplicates


    def print_results(self, results, LaTex = False):
        """
        Print results in either unicode or LaTex.
        :param results: verified results.
        :param LaTex: LaTex printing flag.
        """
        for res_num,res in enumerate(results):
            var_sym = res[0]
            var_gen = lambdify((), var_sym, modules="mpmath")
            a_ = create_series_from_shift_reg(res[3], res[2][:len(res[3])-1], self.depth)
            b_ = (res[1] * (self.depth // len(res[1])))[:self.depth]
            gcf = GeneralizedContinuedFraction(a_,b_)
            if not LaTex:
                print(str(res_num))
                print('lhs: ')
                sympy.pprint(var_sym)
                print('rhs :')
                gcf.print(5)
                print('lhs value: ' + mpmath.nstr(var_gen(), 50))
                print('rhs value: ' + mpmath.nstr(gcf.evaluate(), 50))
            else:
                equation = sympy.Eq(var_sym, gcf.sym_expression(5))
                print(str(res_num)+'. $$ ' + sympy.latex(equation) + ' $$')
                print("\n\n")


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
        self.enum_dps = 50  # working decimal precision for first enumeration
        self.verify_dps = 2000  # working decimal precision for validating results.
        self.lhs_limit = lhs_search_limit
        self.const_sym = sym_constant
        self.const_val = lambdify((), sym_constant, modules="mpmath")
        self.create_an_series = create_series_from_compact_poly
        self.create_bn_series = create_series_from_compact_poly
        if saved_hash == '':
            print('no previous hash table given, initializing hash table...')
            with mpmath.workdps(self.enum_dps):
                start = time()
                self.hash_table = LHSHashTable(
                    range(self.lhs_limit + 1),  # a,b range (allow only non-negative)
                    range(-self.lhs_limit, self.lhs_limit + 1),  # c,d range
                    self.const_val(),  # constant
                    self.threshold)  # length of key
                end = time()
                print('that took {}s'.format(end - start))
        else:
            self.hash_table = LHSHashTable.load_from(saved_hash)

    @staticmethod
    def __number_of_elements(permutation_options: List[List]):
        res = 1
        for l in permutation_options:
            res *= len(l)
        return res

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
        a_coef_list = list(itertools.product(*poly_a))  # all coefficients possibilities for 'a_n'
        neg_poly_b = [[-i for i in b] for b in poly_b]  # for b_n include negative terms
        b_coef_iter = itertools.chain(itertools.product(*poly_b), itertools.product(*neg_poly_b))
        num_iterations = 2 * self.__number_of_elements(poly_b) * self.__number_of_elements(poly_a)
        # create a_n and b_n series fro coefficients.
        an_list = [self.create_an_series(a_coef_list[i], 32) for i in range(len(a_coef_list))]
        # filter out all options resulting in '0' in any series term.
        an_filter = [0 not in an for an in an_list]
        an_list = list(itertools.compress(an_list, an_filter))
        a_coef_list = list(itertools.compress(a_coef_list, an_filter))
        if print_results:
            print('created final enumerations filters after {} s'.format(time() - start))

        counter = 0  # number of permutations passed
        results = []  # list of intermediate results
        start = time()
        for b_coef in b_coef_iter:
            bn = self.create_bn_series(b_coef, 32)
            if 0 in bn:
                continue
            for an_coef in zip(an_list, a_coef_list):
                gcf = EfficientGCF(an_coef[0], bn)  # create gcf from a_n and b_n
                key = int(gcf.evaluate() / self.threshold)  # calculate hash key of gcf value
                if key in self.hash_table:  # find hits in hash table
                    results.append(Match(self.hash_table[key], an_coef[1], b_coef))
                if print_results:
                    counter += 1
                    if counter % 100000 == 0:  # print status.
                        print('passed {} out of {}. found so far {} results'.format(counter, num_iterations,
                                                                                    len(results)))
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
                if mpmath.almosteq(val, t(1), 1 / (self.verify_dps // 20)):
                    # don't keep results that are independent of the constant
                    continue
            except ZeroDivisionError:
                continue

            # create a_n, b_n with huge length, calculate gcf, and verify result.
            an = self.create_an_series(r.rhs_an_poly, 1000)
            bn = self.create_bn_series(r.rhs_bn_poly, 1000)
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
            an = self.create_an_series(r.rhs_an_poly, 1000)
            bn = self.create_bn_series(r.rhs_bn_poly, 1000)
            print_length = max(max(len(r.rhs_an_poly), len(r.rhs_bn_poly)), 5)
            gcf = GeneralizedContinuedFraction(an, bn)
            t = MobiusTransform(r.lhs_coefs)
            sym_lhs = sympy.simplify(t.sym_expression(self.const_sym))
            #print('lhs: {}, an_poly: {}, bn_poly: {}'.format(sym_lhs, r.rhs_an_poly, r.rhs_bn_poly))
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
            # step (2)
            results = self.__first_enumeration(poly_a, poly_b, print_results)
            end = time()
            if print_results:
                print('that took {}s'.format(end - start))
        with mpmath.workdps(self.verify_dps):
            if print_results:
                print('starting to verify results...')
            start = time()
            refined_results = self.__refine_results(results, print_results)  # step (3)
            end = time()
            if print_results:
                print('that took {}s'.format(end - start))
            if print_results:
                self.print_results(refined_results)
        return refined_results


def multi_core_enumeration(sym_constant, lhs_search_limit, saved_hash, poly_a, poly_b, num_cores, splits_size,
                           create_an_series=None, create_bn_series=None, index=0):
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
        if index == (num_cores - 1):
            poly_a[s] = poly_a[s][index * splits_size[s]:]
        else:
            poly_a[s] = poly_a[s][index * splits_size[s]:(index + 1) * splits_size[s]]
    enumerator = EnumerateOverGCF(sym_constant, lhs_search_limit, saved_hash)

    if create_an_series is not None:
        enumerator.create_an_series = create_an_series
    if create_bn_series is not None:
        enumerator.create_bn_series = create_bn_series

    results = enumerator.find_hits(poly_a, poly_b, index == 0)
    enumerator.print_results(results, True)
    return results


def multi_core_enumeration_wrapper(sym_constant, lhs_search_limit, poly_a, poly_b, num_cores, manual_splits_size=None,
                                   saved_hash=None, create_an_series=None, create_bn_series=None):
    """
    a wrapper for enumerating using multi proccessing.
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
    :param create_an_series: a custom function for creating a_n series with poly_a coefficients
    (default is create_series_from_compact_poly)
    :param create_bn_series: a custom function for creating b_n series with poly_b coefficients
    (default is create_series_from_compact_poly)
    :return: results.
    """
    if (saved_hash is None) or (not os.path.isfile(saved_hash)):
        if saved_hash is None:  # if no hash table given, build it here.
            saved_hash = 'tmp_hash.p'
        enumerator = EnumerateOverGCF(sym_constant, lhs_search_limit)
        enumerator.hash_table.save(saved_hash)  # and save it to file (and global instance)
    else:
        if os.name != 'nt':  # if creation of process uses 'Copy On Write' we can benefit from it by
            # loading the hash table to memory here.
            EnumerateOverGCF(sym_constant, lhs_search_limit, saved_hash)

    if manual_splits_size is None:  # naive work split
        manual_splits_size = [len(poly_a[0]) // num_cores]

    # built function for processes
    func = partial(multi_core_enumeration, sym_constant, lhs_search_limit, saved_hash, poly_a, poly_b, num_cores,
                   manual_splits_size, create_an_series, create_bn_series)

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
    final_results = multi_core_enumeration_wrapper(sympy.zeta(2),  # constant to run on
                                                   20,  # lhs limit
                                                   [[i for i in range(12)]] * 3,  # a_n polynomial coefficients
                                                   [[i for i in range(10)]] * 2,  # b_n polynomial coefficients
                                                   2,  # number of cores to run on
                                                   None,  # use naive tiling
                                                   os.path.join('hash_tables', 'zeta2_20_hash.p'),  # if existing
                                                   create_an_series=None,  # use default
                                                   create_bn_series=partial(create_zeta_bn_series, 4)
                                                   )

    # with open('results_of_e_30', 'wb') as file:
    #    pickle.dump(final_results, file)
