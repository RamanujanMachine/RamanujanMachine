import unittest
from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator
from ramanujan.poly_domains.CartesianProductPolyDomain import CartesianProductPolyDomain
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.constants import g_const_dict
from ramanujan.multiprocess_enumeration import multiprocess_enumeration


def get_testable_data(refined_res_list):
    """
    On those test we have no use for the bloom key or lhs dict pointers
    extract the meaning full data for testing
    """
    data = []
    for i in refined_res_list:
        data.append((
            i.rhs_an_poly,
            i.rhs_bn_poly,
            i.c_top,
            i.c_bot))
    return data


class APITests(unittest.TestCase):
    def test_MITM_api1(self):
        lhs = LHSHashTable('e_lhs_dept5_db', 5, [g_const_dict['e']])

        poly_search_domain = CartesianProductPolyDomain(
            1, [-5, 5],
            1, [-5, 5])

        enumerator = EfficientGCFEnumerator(lhs, poly_search_domain, [g_const_dict['e']])

        results = get_testable_data(enumerator.full_execution())

        self.assertEqual(len(results), 20)
        self.assertIn(
            ((4, 2), (0, 1), (1, 1), (-1, 1)),
            results)
        self.assertIn(
            ((1, 1), (1, 0), (1, 0), (-2, 1)),
            results)

    def test_MITM_api2(self):
        lhs = LHSHashTable('zeta3.lhs.dept14.db', 14, [g_const_dict['zeta'](3)])

        poly_search_domain = Zeta3Domain1(
            [(2, 2), (1, 1), (1, 17), (1, 5)],  # an coefs
            (-16, -1)  # bn coef
            )

        enumerator = EfficientGCFEnumerator(
            lhs,
            poly_search_domain,
            [g_const_dict['zeta'](3)]
            )

        results = get_testable_data(enumerator.full_execution())

        self.assertEqual(len(results), 3)
        self.assertIn(
            ((2, 1, 3, 1), (-1,), (8, 0), (0, 7)),
            results)
        self.assertIn(
            ((2, 1, 5, 2), (-16,), (12, 0), (0, 7)),
            results)
        self.assertIn(
            ((2, 1, 17, 5), (-1,), (6, 0), (0, 1)),
            results)

    def test_MITM_api3(self):
        saved_hash = 'pi_lhs_dept20'
        lhs_search_limit = 20
        lhs = LHSHashTable(
            saved_hash,
            lhs_search_limit,
            [g_const_dict['pi']])

        poly_search_domain = CartesianProductPolyDomain(
            1, [-13, 13],
            2, [-11, 11])

        # create an enumerator to iter thought the poly domain and compare it to the
        # lhs table
        enumerator = EfficientGCFEnumerator(
            lhs,
            poly_search_domain,
            [g_const_dict['pi']]
            )

        results = enumerator.full_execution()
        self.assertEqual(len(results), 46)

    def test_MITM_multiprocessing(self):
        """
        This is the same test as test_MITM_api2, but using multiprocessing and expanding the search domain.
        This causes some 'code duplicate' warnings.
        This tested domain is large, and will scan duplicates (two scaled results). This is not a mistake, since we're
        interested in the way the poly domain is split, and we want to test that no part of the domain is missed out.
        """
        lhs = LHSHashTable('zeta3_lhs_dept20.db', 20, [g_const_dict['zeta'](3)])

        poly_search_domain = Zeta3Domain1(
            [(2, 2), (1, 1), (1, 100), (1, 100)],  # an coefs
            (-50, -1)  # bn coef
            )

        results = multiprocess_enumeration(
            EfficientGCFEnumerator,
            lhs,
            poly_search_domain,
            [g_const_dict['zeta'](3)],
            4)

        results = get_testable_data(results)

        self.assertEqual(len(results), 7)
        self.assertIn(
            ((2, 1, 3, 1), (-1,), (8, 0), (0, 7)),
            results)
        self.assertIn(
            ((2, 1, 5, 2), (-16,), (12, 0), (0, 7)),
            results)
        self.assertIn(
            ((2, 1, 6, 2), (-4,), (16, 0), (0, 7)),
            results)
        self.assertIn(
            ((2, 1, 17, 5), (-1,), (6, 0), (0, 1)),
            results)
        self.assertIn(
            ((2, 1, 21, 7), (-49,), (8, 0), (0, 1)),
            results)
        self.assertIn(
            ((2, 1, 34, 10), (-4,), (12, 0), (0, 1)),
            results)
        self.assertIn(
            ((2, 1, 51, 15), (-9,), (18, 0), (0, 1)),
            results)

    def test_poly_domain_split(self):
        """
        making sure that the domain is splitted currectly
        1. checking if the approximate size of a splitted domain is the same as the
           original domain
        2. checking if all the items in the splitted domain are
        """
        # coef ranges are 29, aiming for primes that screw with even splitting
        cartesian_domain = CartesianProductPolyDomain(
            2, [-30, 30],
            3, [-10, 19])
        splitted_cartesian_domain = cartesian_domain.split_domains_to_processes(5)

        self.assertEqual(cartesian_domain.num_iterations, sum([i.num_iterations for i in splitted_cartesian_domain]))
        print('passed1')
        # the zeta domain checks for convergences condition when iterating over coefs
        # so the approximate size is bigger then the one actually used
        original_zeta_domain = Zeta3Domain1(
            [(2, 10), (1, 1), (1, 30), (1, 10)],
            (-10, -1)
        )
        splitted_zeta_domain = original_zeta_domain.split_domains_to_processes(7)

        all_zeta_polys = [i for i in original_zeta_domain.iter_polys('b')]
        for sub_domain in splitted_zeta_domain:
            for polys in sub_domain.iter_polys('a'):
                # checking if a value is present this way imporves execution times drasticly
                try:
                    all_zeta_polys.remove(polys)
                except ValueError:
                    self.assertIn(polys, all_zeta_polys)

        self.assertEqual(len(all_zeta_polys), 0)
if __name__ == '__main__':
    unittest.main()
