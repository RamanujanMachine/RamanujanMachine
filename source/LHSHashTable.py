import sympy
import struct
import mpmath
import pickle
import itertools
from sympy import lambdify
from datetime import datetime
from time import time
from pybloom_live import BloomFilter
from functools import partial, reduce
from math import gcd

class LHSHashTable(object):
    @staticmethod
    def are_co_prime(integers):
        common = integers[-1]
        for x in integers:
            common = gcd(x, common)
            if common == 1:
                return True
        return False

    @staticmethod
    def prod(coefs, consts):
        ret = coefs[0]
        for i in range(len(coefs) - 1):
            ret += consts[i] * coefs[i + 1]
        return ret

    @staticmethod
    def lhs_hash_name_to_shelve_name(name):
        return name.split('.')[0] + '.db'

    def _get_by_key(self, key):
        if self.s is None:
            with open(self.s_name, 'rb') as f:
                self.s = pickle.load(f)
        vals = struct.unpack(self.pack_format, self.s[str(key)])
        return vals[:self.n_constants], vals[-self.n_constants:]

    def evaluate(self, key, constant_values):
        c_top, c_bottom = self._get_by_key(key)
        numerator = self.prod(c_top, constant_values)
        denominator = self.prod(c_bottom, constant_values)
        return mpmath.mpf(numerator) / mpmath.mpf(denominator)

    def evaluate_sym(self, key, symbols):
        c_top, c_bottom = self._get_by_key(key)
        numerator = self.prod(c_top, symbols)
        denominator = self.prod(c_bottom, symbols)
        return numerator / denominator

    def __init__(self, name, search_range, const_vals, threshold) -> None:
        """
        hash table for LHS. storing values in the form of (a + b*x_1 + c*x_2 + ...)/(d + e*x_1 + f*x_2 + ...)
        :param search_range: range for value coefficient values
        :param const_vals: constants for x.
        :param threshold: decimal threshold for comparison. in fact, the keys for hashing will be the first
                            -log_{10}(threshold) digits of the value. for example, if threshold is 1e-10 - then the
                            first 10 digits will be used as the hash key.
        """
        self.name = name
        self.s_name = self.lhs_hash_name_to_shelve_name(name)
        self.threshold = threshold
        key_factor = 1 / threshold
        self.max_key_length = len(str(int(key_factor))) * 2

        # create blacklist of rational numbers
        coef_possibilities = [i for i in range(-search_range, search_range + 1)]
        coef_possibilities.remove(0)
        rational_options = itertools.product(*[coef_possibilities, coef_possibilities])
        rational_keys = [int((mpmath.mpf(ratio[0]) / ratio[1]) * key_factor) for ratio in rational_options]
        # +-1 for numeric errors in keys.
        rational_blacklist = set(rational_keys + [x + 1 for x in rational_keys] + [x - 1 for x in rational_keys])

        # create enumeration lists
        constants = [mpmath.mpf(1)] + const_vals
        self.n_constants = len(constants)
        coefs_top = [range(-search_range, search_range + 1)] * len(constants)  # numerator range
        coefs_bottom = [range(-search_range, search_range + 1)] * len(constants)  # denominator range
        coef_top_list = itertools.product(*coefs_top)
        coef_bottom_list = list(itertools.product(*coefs_bottom))
        denominator_list = [sum(i * j for (i, j) in zip(c_bottom, constants)) for c_bottom in coef_bottom_list]

        # start enumerating
        t = time()

        self.max_capacity = (search_range * 2 + 1) ** (len(constants) * 2)
        self.s = {}
        self.pack_format = 'll' * self.n_constants
        self.bloom = BloomFilter(capacity=self.max_capacity, error_rate=0.05)
        for c_top in coef_top_list:
            numerator = sum(i * j for (i, j) in zip(c_top, constants))
            if numerator <= 0:  # allow only positive values to avoid duplication
                continue
            numerator = mpmath.mpf(numerator)
            for c_bottom, denominator in zip(coef_bottom_list, denominator_list):
                if reduce(gcd, c_top + c_bottom) != 1:  # avoid expressions that can be simplified easily
                    continue
                if denominator == 0:  # don't store inf or nan.
                    continue
                val = numerator / denominator
                key = int(val * key_factor)
                if key in rational_blacklist:
                    # don't store values that are independent of the constant (e.g. rational numbers)
                    continue
                str_key = str(key)
                self.s[str_key] = struct.pack(self.pack_format, *[*c_top, *c_bottom])  # store key and transformation
                self.bloom.add(str_key)

        with open(self.s_name, 'wb') as f:
            pickle.dump(self.s, f)
        self.s = None
        print('initializing LHS dict: {}'.format(time() - t))

    def __contains__(self, item):
        """
        operator 'in'
        :param item: key
        :return: true of false
        """
        return item in self.bloom

    def __getitem__(self, item):
        """
        operator []
        :param item: key
        :return: transformation of x
        """
        return self._get_by_key(item)

    def __eq__(self, other):
        """
        operator ==
        :param other: other hash table.
        :return:
        """
        if type(other) != type(self):
            return False
        ret = self.threshold == other.threshold
        # ret &= sorted(self.s.keys()) == sorted(other.s.keys())
        return ret

    def save(self):
        """
        save the hash table as file
        """
        with open(self.name, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load_from(cls, name):
        """
        load hash table from file (or global instance)
        :param name:
        :return:
        """
        with open(name, 'rb') as f:
            ret = pickle.load(f)
        ret.s_name = ret.lhs_hash_name_to_shelve_name(name)
        return ret