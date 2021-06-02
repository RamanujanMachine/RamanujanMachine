from .AbstractPolyDomains import AbstractPolyDomains
from ..utils.utils import iter_series_items_from_compact_poly
from itertools import product
from copy import deepcopy
from numpy import array_split


class CartesianProductPolyDomain(AbstractPolyDomains):
    """
    This poly domain will generate all combinations for a(n) and b(n) coefs without complex dependence between the two
    """
    def __init__(self, a_deg, a_coef_range, b_deg, b_coef_range, an_leading_coef_positive=True, *args, **kwargs):
        """
        If all of an's coefs can get both positive and negative values, then we might get two iterations for any set of
        coefs, with opposite signs. Those two series will converge to the same value, but with a different sign, hence
        it is a redundant run we can skip. Based on an_leading_coef_positive we will try to detect those cases and skip
        them
        """
        self.a_deg = a_deg
        # expanding the range to a different range for each coef
        # allows us to use the same functions for decedent classes
        self.a_coef_range = [list(a_coef_range) for _ in range(a_deg + 1)]
        if an_leading_coef_positive and self.a_coef_range[-1][0] <= 0:
            self.a_coef_range[0][0] = 1

        self.b_deg = b_deg
        self.b_coef_range = [b_coef_range for _ in range(b_deg + 1)]

        self._setup_metadata()
        super().__init__()

    def _setup_metadata(self):
        """
        This function generates and stores values that should not change throughout the run.
        It continues __init__'s job, but holds code that is used in classes that extend this class, so it was
        moved to a separate function.
        """
        self.an_length = self.get_an_length()
        self.bn_length = self.get_bn_length()
        self.num_iterations = self.an_length * self.bn_length

        self.an_domain_range, self.bn_domain_range = self.dump_domain_ranges()

    @staticmethod
    def _range_size(coef_range):
        return coef_range[1] - coef_range[0] + 1

    @staticmethod
    def domain_size_by_var_ranges(var_ranges):
        size = 1
        for var_range in var_ranges:
            size *= CartesianProductPolyDomain._range_size(var_range)
        return size

    @staticmethod
    def expand_coef_range_to_full_domain(coef_ranges):
        return [[i for i in range(coef[0], coef[1] + 1)] for coef in coef_ranges]

    def get_an_length(self):
        return CartesianProductPolyDomain.domain_size_by_var_ranges(self.a_coef_range)

    def get_bn_length(self):
        return CartesianProductPolyDomain.domain_size_by_var_ranges(self.b_coef_range)

    @staticmethod
    def get_calculation_method():
        # both an and bn are regular compact polys
        return iter_series_items_from_compact_poly, \
            iter_series_items_from_compact_poly

    def dump_domain_ranges(self):
        an_domain = CartesianProductPolyDomain.expand_coef_range_to_full_domain(self.a_coef_range)
        bn_domain = CartesianProductPolyDomain.expand_coef_range_to_full_domain(self.b_coef_range)

        return an_domain, bn_domain

    def iter_polys(self, primary_looped_domain):
        an_domain, bn_domain = self.dump_domain_ranges()

        if primary_looped_domain == 'a':
            a_coef_iter = product(*an_domain)
            for a_coef in a_coef_iter:
                b_coef_iter = product(*bn_domain)
                for b_coef in b_coef_iter:
                    yield a_coef, b_coef
        else:
            b_coef_iter = product(*bn_domain)
            for b_coef in b_coef_iter:
                a_coef_iter = product(*an_domain)
                for a_coef in a_coef_iter:
                    yield a_coef, b_coef

    def get_a_coef_iterator(self):
        return product(*self.an_domain_range)

    def get_b_coef_iterator(self):
        return product(*self.bn_domain_range)

    def get_individual_polys_generators(self):
        # for backwards compatibility.
        return self.get_a_coef_iterator(), self.get_b_coef_iterator()

    @staticmethod
    def _get_metadata_on_var_ranges(ranges, series):
        ranges_metadata = []
        for i, v in enumerate(ranges):
            ranges_metadata.append({
                'range': v,
                'size': v[1]-v[0]+1,
                'series': series,
                'index': i
                })
        return ranges_metadata

    def split_domains_to_processes(self, number_of_instances):
        """
        When using multiprocessing, we'll split the search domain to several polyDomain and iter over each one
        in a different process.
        This function will split the domain to number_of_instances sub-domains. To do so, we'll find the coef
        with the biggest range, and split it as evenly as possible to different instances.
        """
        all_coef_ranges = self._get_metadata_on_var_ranges(self.a_coef_range, 'a')
        all_coef_ranges += self._get_metadata_on_var_ranges(self.b_coef_range, 'b')

        biggest_range = max(all_coef_ranges, key=lambda x: x['size'])

        # when splitting over a huge number of processes (probably over different clients) we'll not be able to
        # split the domain using only one coef. If that's the case, we'll just split the biggest coef as much as 
        # we can, and use the same logic recursively over every sub-domain.
        number_of_sub_arrays = min(number_of_instances, biggest_range['size'])
        
        sub_domains = []
        for chunk_items in array_split(range(biggest_range['range'][0], biggest_range['range'][1] + 1),
                                       number_of_sub_arrays):
            chunk_range = [int(chunk_items[0]), int(chunk_items[-1])]
            next_instance = deepcopy(self)
            if biggest_range['series'] == 'a':
                next_instance.a_coef_range[biggest_range['index']] = chunk_range
            else:
                next_instance.b_coef_range[biggest_range['index']] = chunk_range

            next_instance._setup_metadata()
            sub_domains.append(next_instance)

        if biggest_range['size'] < number_of_instances:
            # divide the required instances over all of the sub domains
            smaller_sub_domains = []
            for i, sub_domain in zip(array_split(range(number_of_instances), len(sub_domains)), sub_domains):
                smaller_sub_domains += sub_domain.split_domains_to_processes(len(i))

            return smaller_sub_domains

        return sub_domains
