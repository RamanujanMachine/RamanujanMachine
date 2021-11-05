from .AbstractPolyDomains import AbstractPolyDomains
from ..utils.utils import iter_series_items_from_compact_poly
from itertools import product
from copy import deepcopy
from numpy import array_split

CHECKPOINT_DUMP_SIZE = 5_000
ALLOW_LOWER_DEGREE = False


class CartesianProductPolyDomain(AbstractPolyDomains):
    """
    This poly domain will generate all combinations for a(n) and b(n) coefficients without complex dependence between
    the two
    """
    def __init__(self, a_deg, a_coef_range, b_deg, b_coef_range, an_leading_coef_positive=True,
                 only_balanced_degrees=False, use_strict_convergence_cond=False, *args, **kwargs):
        """
        a_deg - an's polynomial degree
        a_coef_range - The coefficient range iterated for every coefficient in an
        b_deg - bn's polynomial degree
        b_coef_range - The coefficient range iterated for every coefficient in bn
        an_leading_coef_positive - A GCF can be inflated by multiplying an*c and bn*c^2. The inflation will converge to
            the same value (up to multiplying by a factor). By forcing an_leading coefficient to be positive we can
            remove cases when c<0 and make identifying those cases more easy
        only_balanced_degrees - forces deg(an)*2=deg(bn). Read ramanujan machine paper for more information about
            convergence conditions.
        use_strict_convergence_cond - discard cases when discriminate = 0. Read ramanujan machine paper for more
            information about convergence conditions.
        """
        self.a_deg = a_deg
        # expanding the range to a different range for each coefficient
        # allows us to use the same functions for decedent classes
        self.a_coef_range = [list(a_coef_range) for _ in range(a_deg + 1)]
        if an_leading_coef_positive and self.a_coef_range[-1][0] <= 0:
            self.a_coef_range[0][0] = 1

        self.b_deg = b_deg
        self.b_coef_range = [b_coef_range for _ in range(b_deg + 1)]
        self.only_balanced_degrees = only_balanced_degrees
        self.use_strict_convergence_cond = use_strict_convergence_cond

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

    def expand_coef_range_to_full_domain(self, coef_ranges):
        domain = [[i for i in range(coef[0], coef[1] + 1)] for coef in coef_ranges]
        if self.only_balanced_degrees and 0 in domain[0]:
            domain[0].remove(0)
        return domain

    def get_an_length(self):
        return self.domain_size_by_var_ranges(self.a_coef_range)

    def get_bn_length(self):
        return self.domain_size_by_var_ranges(self.b_coef_range)

    def get_an_degree(self, an_coefs):
        if ALLOW_LOWER_DEGREE:
            deg = self.a_deg
            for i in an_coefs:
                if i == 0:
                    deg -= 1
                else:
                    break
            return deg

        return self.a_deg

    def get_bn_degree(self, bn_coefs):
        if ALLOW_LOWER_DEGREE:
            deg = self.b_deg
            for i in bn_coefs:
                if i == 0:
                    deg -= 1
                else:
                    break
            return deg

        return self.b_deg

    @staticmethod
    def _get_compact_poly_deg(coefs):
        deg = len(coefs) - 1
        for coef in coefs:
            if coef != 0:
                return deg
            else:
                deg -= 1

    @staticmethod
    def get_calculation_method():
        # both an and bn are regular compact polys
        return iter_series_items_from_compact_poly, iter_series_items_from_compact_poly

    def dump_domain_ranges(self):
        an_domain = self.expand_coef_range_to_full_domain(self.a_coef_range)
        bn_domain = self.expand_coef_range_to_full_domain(self.b_coef_range)

        return an_domain, bn_domain

    def filter_gcfs(self, an_coefs, bn_coefs):
        """
        Some GCFs will not converge, and some are duplicates of other GCFs.
        This function filter these cases out.
        """
        # For un-balanced degrees, we have no filtering conditions
        if (len(an_coefs) - 1) * 2 != len(bn_coefs) - 1:
            return self.only_balanced_degrees

        # Discard non-converging cases
        if 4 * bn_coefs[0] < -1 * (an_coefs[0]**2):
            return False

        # On equality in convergence condition, some cases converge and some not 
        if self.use_strict_convergence_cond and 4 * bn_coefs[0] == -1 * (an_coefs[0]**2):
            return False

        return True

    def iter_polys(self, primary_looped_domain):
        """
        This function iterate pairs of an and bn coefficients from the domain.
        Some enumerators cache series items, and primary_looped_domain is used to determine the nested loop order that 
        fit the caching mechanism. Only the nested series needs to be cached.
        The outer looped series is called pn, and the inner series sn.
        """
        def _get_coefs_in_order():
            # Helper function to order pn and sn back to an and bn 
            if primary_looped_domain == 'a':
                return pn_coef, sn_coef
            else:
                return sn_coef, pn_coef

        # setting pn and sn from original series
        if primary_looped_domain == 'a':
            pn_coef_range = self.a_coef_range
            sn_coef_range = self.b_coef_range
        else:
            pn_coef_range = self.b_coef_range
            sn_coef_range = self.a_coef_range

        pn_domain = self.expand_coef_range_to_full_domain(pn_coef_range)
        sn_domain = self.expand_coef_range_to_full_domain(sn_coef_range)

        # Discard items that we're already calculated from the outer series
        pn_iterator = product(*pn_domain)

        items_passed = 0
        for pn_coef in pn_iterator:
            sn_iterator = product(*sn_domain)           
            for sn_coef in sn_iterator:
                if self.filter_gcfs(*_get_coefs_in_order()):
                    yield _get_coefs_in_order()
                items_passed += 1

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
        This function will split the domain to number_of_instances sub-domains. To do so, we'll find the coefficient
        with the biggest range, and split it as evenly as possible to different instances.
        """
        all_coef_ranges = self._get_metadata_on_var_ranges(self.a_coef_range, 'a')
        all_coef_ranges += self._get_metadata_on_var_ranges(self.b_coef_range, 'b')

        biggest_range = max(all_coef_ranges, key=lambda x: x['size'])

        # when splitting over a huge number of processes (probably over different clients) we'll not be able to
        # split the domain using only one coefficient. If that's the case, we'll just split the biggest coefficient as
        # much as we can, and use the same logic recursively over every sub-domain.
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
