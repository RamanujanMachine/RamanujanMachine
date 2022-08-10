import os
import struct
import pickle
import mpmath
import itertools
from time import time
from pybloom_live import BloomFilter
from functools import reduce
from math import gcd
from ramanujan.utils.utils import create_mpf_const_generator

from ramanujan.constants import g_N_initial_search_dps

# precision required from table
DEFAULT_THRESHOLD = 10**-10


class LHSHashTable(object):
    """
    This class makes use of bloom filters and a regular representation to improve performance
    LHS items are stored in their "raw" form on a file called self.s_name. This file is only
    opened when needed, to reduce memory consumptions.
    The bloom filter is always loaded and used to determine if a LHS value is in the database
    all LHS possibilities within computed domain
    """

    def __init__(
        self, name, search_range, const_vals, threshold=DEFAULT_THRESHOLD
    ) -> None:
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
        self.constant_generator = create_mpf_const_generator(const_vals)
        const_vals = [const() for const in self.constant_generator]
        constants = [mpmath.mpf(1)] + const_vals
        self.n_constants = len(constants)

        self.max_capacity = (search_range * 2 + 1) ** (self.n_constants * 2)
        self.pack_format = "ll" * self.n_constants
        self.lhs_possibilities = {}
        self.bloom = BloomFilter(capacity=self.max_capacity, error_rate=0.05)

        start_time = time()

        if os.path.isfile(self.s_name):
            print(f"loading from {self.s_name}")
            self._load_from_file(self.s_name)
        else:
            print("no existing db found, generating dict")
            with mpmath.workdps(g_N_initial_search_dps):
                self._enumerate_lhs_domain(constants, search_range, key_factor)

        with open(self.s_name, "wb") as f:
            pickle.dump(self.lhs_possibilities, f)

        # after init, deleteing self.lhs_possibilities to free unused memory
        self.lhs_possibilities = None
        print("initializing LHS dict: {}".format(time() - start_time))

    @staticmethod
    def _create_rational_numbers_blacklist(search_range, key_factor):
        # LHS numerator and denominator might cancel out and LHS will be rational.
        # Those options are not relevant
        coef_possibilities = [i for i in range(-search_range, search_range + 1)]
        coef_possibilities.remove(0)
        rational_options = itertools.product(coef_possibilities, coef_possibilities)
        rational_keys = [
            int((mpmath.mpf(num) / denom) * key_factor)
            for num, denom in rational_options
        ]
        # +-1 for numeric errors in keys.
        return set(
            rational_keys
            + [x + 1 for x in rational_keys]
            + [x - 1 for x in rational_keys]
        )

    def _load_from_file(self, db_path):
        with open(db_path, "rb") as f:
            self.lhs_possibilities = pickle.load(f)
        for key in self.lhs_possibilities.keys():
            self.bloom.add(key)

    def _enumerate_lhs_domain(self, constants, search_range, key_factor):
        rational_blacklist = LHSHashTable._create_rational_numbers_blacklist(
            search_range, key_factor
        )

        # Create enumeration lists
        coefs_top = [
            range(-search_range, search_range + 1)
        ] * self.n_constants  # numerator range
        coefs_bottom = [
            range(-search_range, search_range + 1)
        ] * self.n_constants  # denominator range
        coef_top_list = itertools.product(*coefs_top)
        coef_bottom_list = list(itertools.product(*coefs_bottom))
        denominator_list = [
            sum(i * j for (i, j) in zip(c_bottom, constants))
            for c_bottom in coef_bottom_list
        ]

        # start enumerating
        for c_top in coef_top_list:
            numerator = sum(i * j for (i, j) in zip(c_top, constants))
            if numerator <= 0:  # allow only positive values to avoid duplication
                continue
            numerator = mpmath.mpf(numerator)

            for c_bottom, denominator in zip(coef_bottom_list, denominator_list):
                if (
                    reduce(gcd, c_top + c_bottom) != 1
                ):  # avoid expressions that can be simplified easily
                    continue
                if denominator == 0:  # don't store inf or nan.
                    continue
                val = numerator / denominator
                key = int(val * key_factor)

                if key in rational_blacklist:
                    # don't store values that are independent of the constant (e.g. rational numbers)
                    continue

                str_key = str(key)
                self._add_to_lhs_possibilities(str_key, c_top, c_bottom)
                self.bloom.add(str_key)

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
        return ret

    @staticmethod
    def lhs_hash_name_to_shelve_name(name):
        return name.split(".")[0] + ".db"

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

    def _get_by_key(self, key):
        with open(self.s_name, "rb") as f:
            if self.lhs_possibilities is None:
                self.lhs_possibilities = pickle.load(f)
            values = []
            for lhs_match in self.lhs_possibilities[str(key)]:
                vals = struct.unpack(self.pack_format, lhs_match)
                values.append([vals[: self.n_constants], vals[-self.n_constants :]])
            return values

    @classmethod
    def load_from(cls, name):
        """
        load hash table from file (or global instance)
        :param name:
        :return:
        """
        with open(name, "rb") as f:
            ret = pickle.load(f)
        ret.s_name = ret.lhs_hash_name_to_shelve_name(name)
        return ret

    def _add_to_lhs_possibilities(self, str_key, c_top, c_bottom):
        # store key and transformation
        if str_key in self.lhs_possibilities:
            self.lhs_possibilities[str_key].append(
                struct.pack(self.pack_format, *[*c_top, *c_bottom])
            )
        else:
            self.lhs_possibilities[str_key] = [
                struct.pack(self.pack_format, *[*c_top, *c_bottom])
            ]

    def evaluate(self, key):
        # this function will usually be called under a different mpf workdps
        # generating constant_values again here to match scope's precision
        const_vals = [const() for const in self.constant_generator]
        stored_values = self._get_by_key(key)
        evaluated_values = []
        for c_top, c_bottom in stored_values:
            numerator = self.prod(c_top, const_vals)
            denominator = self.prod(c_bottom, const_vals)
            evaluated_values.append(
                (mpmath.mpf(numerator) / mpmath.mpf(denominator), c_top, c_bottom)
            )

        return evaluated_values

    def evaluate_sym(self, key, symbols):
        stored_values = self._get_by_key(key)
        evaluated_values = []

        for c_top, c_bottom in stored_values:
            numerator = self.prod(c_top, symbols)
            denominator = self.prod(c_bottom, symbols)
            evaluated_values.append(numerator / denominator)
        return evaluated_values

    def save(self):
        """
        save the hash table as file
        """
        with open(self.name, "wb") as f:
            pickle.dump(self, f)
