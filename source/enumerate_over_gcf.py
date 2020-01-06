import pickle
import mpmath
import numpy as np
from sympy import lambdify
import sympy
from time import time
import itertools
from massey import create_series_from_compact_poly
from mobius import GeneralizedContinuedFraction, MobiusTransform
import sys
from collections import namedtuple
from typing import List

Match = namedtuple('Match', 'lhs_coefs rhs_an_poly rhs_bn_poly')


class LHSHashTable(object):
    def __init__(self, limit, const_val, threshold) -> None:
        self.s = {}
        self.threshold = threshold
        for a in range(-limit, limit):
            for b in range(-limit, limit):
                for c in range(-limit, limit):
                    for d in range(-limit, limit):
                        denom = c*const_val + d
                        if denom == 0:
                            continue
                        val = (a*const_val + b) / denom
                        if mpmath.isnan(val) or mpmath.isinf(val):
                            continue
                        key = int(mpmath.nint(val / self.threshold))
                        if key in self.s:
                            continue
                        self.s[key] = np.array([[a, b], [c, d]], dtype=object)

    def __contains__(self, item):
        return item in self.s

    def __getitem__(self, item):
        return self.s[item]

    def __eq__(self, other):
        ret = type(other) == type(self)
        ret &= self.threshold == other.threshold
        ret &= sorted(self.s.keys()) == sorted(other.s.keys())
        return ret

    def save(self, name):
        with open(name, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load_from(cls, name):
        with open(name, 'rb') as f:
            ret = pickle.load(f)
        return ret


class EnumerateOverGCF(object):
    def __init__(self, sym_constant, lhs_search_limit, rhs_search_limit, saved_hash=''):
        self.threshold = 1e-10
        self.enum_dps = 50
        self.verify_dps = 200
        self.lhs_limit = lhs_search_limit
        self.rhs_limit = rhs_search_limit
        self.const_sym = sym_constant
        self.const_val = lambdify((), sym_constant, modules="mpmath")
        if saved_hash == '':
            print('no previous hash table given, initializing hash table...')
            with mpmath.workdps(self.enum_dps):
                start = time()
                self.hash_table = LHSHashTable(self.lhs_limit, self.const_val(), self.threshold)
                end = time()
                print('that took {}s'.format(end-start))
        else:
            self.hash_table = LHSHashTable.load_from(saved_hash)

    def create_poly_enumeration(self, poly_deg_a, poly_deg_b):
        coefs = [[i for i in range(-self.rhs_limit, self.rhs_limit+1)]]
        a = coefs*poly_deg_a
        b = coefs*poly_deg_b
        a_polynomials = list(itertools.product(*a))
        b_polynomials = list(itertools.product(*b))
        return a_polynomials, b_polynomials

    def create_SR_enumeration(self, SR_len_a, SR_len_b, reg_size):
        rng = range(-self.rhs_limit, self.rhs_limit+1)
        cart_prod = list(itertools.product(rng, repeat=reg_size))  # list of tuples
        polynoms = [list(prod) for prod in cart_prod]  # lists of lists
        registers = [list(prod) for prod in cart_prod]  # lists of lists
        for poly in polynoms:
            poly.insert(0, 1)
        domain = list(itertools.product(polynoms, polynoms, registers, registers))  # list of tuples
        return domain

    def first_enumeration(self, poly_deg_a, poly_deg_b, print_results=True):
        results = []
        counter = 0
        a_coef_list, b_coef_list = self.create_enumeration(poly_deg_a, poly_deg_b)
        n_iterations = len(a_coef_list) * len(b_coef_list)
        for b_coef in b_coef_list:
            bn = create_series_from_compact_poly(b_coef, 100)
            if any(b == 0 for b in bn):
                continue
            for a_coef in a_coef_list:
                an = create_series_from_compact_poly(a_coef, 100)
                if any(a == 0 for a in an[1:]):
                    continue
                counter += 1
                if (counter % 1000) == 0 and print_results:
                    print('passed {} permutations out of {}. found so far {} matches'.format(
                        counter, n_iterations, len(results)))
                try:
                    gcf = GeneralizedContinuedFraction(an, bn)
                    key = int(mpmath.nint(gcf.evaluate() / self.threshold))
                except ZeroDivisionError:
                    continue
                except:
                    print("Unexpected error:", sys.exc_info()[0])
                    print('happened with:\n an_poly={}, bn_poly={}'.format(a_coef, b_coef))
                    continue
                if key in self.hash_table:
                    results.append(Match(self.hash_table[key], a_coef, b_coef))
        return results

    def refine_results(self, intermediate_results: List[Match], print_results=True):
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
                if mpmath.isinf(val) or mpmath.isnan(val):
                    continue
                if mpmath.almosteq(val, t(1), 1/(self.verify_dps//20)):
                    continue
            except ZeroDivisionError:
                continue
            an = create_series_from_compact_poly(r.rhs_an_poly, 1000)
            bn = create_series_from_compact_poly(r.rhs_bn_poly, 1000)
            gcf = GeneralizedContinuedFraction(an, bn)
            val_str = mpmath.nstr(val, 100)
            rhs_str = mpmath.nstr(gcf.evaluate(), 100)
            if val_str == rhs_str:
                results.append(r)
        return results

    def print_results(self, results: List[Match]):
        for r in results:
            an = create_series_from_compact_poly(r.rhs_an_poly, 1000)
            bn = create_series_from_compact_poly(r.rhs_bn_poly, 1000)
            gcf = GeneralizedContinuedFraction(an, bn)
            t = MobiusTransform(r.lhs_coefs)
            print('lhs: ')
            sym_lhs = t.sym_expression(self.const_sym)
            sympy.pprint(sympy.simplify(sym_lhs))
            print('rhs: ')
            gcf.print(5)
            print('lhs value: ' + mpmath.nstr(t(self.const_val()), 50))
            print('rhs value: ' + mpmath.nstr(gcf.evaluate(), 50))

    def find_hits(self, poly_deg_a, poly_deg_b, print_results=True):
        with mpmath.workdps(self.enum_dps):
            if print_results:
                print('starting preliminary search...')
            start = time()
            results = self.first_enumeration(poly_deg_a, poly_deg_b, print_results)
            end = time()
            if print_results:
                print('that took {}s'.format(end-start))
        with mpmath.workdps(self.verify_dps):
            if print_results:
                print('starting to verify results...')
            start = time()
            refined_results = self.refine_results(results, print_results)
            end = time()
            if print_results:
                print('that took {}s'.format(end-start))
            if print_results:
                self.print_results(refined_results)
        return refined_results


#enumerator = EnumerateOverGCF(sympy.pi, 20, 3, saved_hash='pi20_hash.p')
#enumerator = EnumerateOverGCF(sympy.pi, 20, 10)
#enumerator.hash_table.save('pi20_hash.p')
#results = enumerator.find_hits(2, 3)
