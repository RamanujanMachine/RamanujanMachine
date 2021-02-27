import unittest
from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator
from ramanujan.poly_domains.CartesianProductPolyDomain import CartesianProductPolyDomain
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.constants import g_const_dict


def get_testable_data(refined_res_list):
    '''
    On those test we have no use for the bloom key or lhs dict pointers
    exctact the meaning full data for testing
    '''
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
            [(2, 2), (1, 1), (1, 17), (1, 5)], # an coefs
            (-16, -1) # bn coef
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
            1, [-13,13],
            2, [-11,11])

        # create an enumerator to iter thought the poly domain and compare it to the
        # lhs table
        enumerator = EfficientGCFEnumerator(
            lhs,
            poly_search_domain,
            [g_const_dict['pi']]
            )

        results = enumerator.full_execution()
        self.assertEqual(len(results), 46)


if __name__ == '__main__':
    unittest.main()
