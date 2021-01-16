import unittest
import sympy
from ramanujan.LHSHashTable import *
from ramanujan.enumerators.EfficentGCFEnumerator import *
from ramanujan.poly_domains.CartesianProductPolyDomain import *
from ramanujan.poly_domains.Zeta3Domain1 import *

def get_testable_data(refind_res_list):
    '''
        On those test we have no use for the bloom key or lhs dict pointers
        exctact the meaning full data for testing
    '''
    data = []
    for i in refind_res_list:
        data.append((
            i.rhs_an_poly,
            i.rhs_bn_poly,
            i.c_top,
            i.c_bot))
    return data


class APITests(unittest.TestCase):

    def test_MITM_api1(self):
        saved_hash = 'e_lhs_dept5_db'
        lhs_search_limit = 5
        lhs = LHSHashTable(saved_hash, lhs_search_limit, [g_const_dict['e']]) 

        poly_search_domain = CartesianProductPolyDomain(
            2, [-5,5], 
            2, [-5,5]) 

        enumerator = EfficentGCFEnumerator(lhs, poly_search_domain, [g_const_dict['e']],
            lhs_search_limit)

        results = get_testable_data(enumerator.full_execution())

        self.assertEqual(len(results), 42)
        self.assertIn(
            ((0, 4, 2), (0, 0, 1), (1, 1), (-1,1)), 
            results)
        self.assertIn(
            ((0, 1, 1), (0, 1, 0), (1, 0), (-2,1)), 
            results)

    def test_MITM_api2(self):
        saved_hash = 'zeta3.lhs.dept14.db'
        lhs_search_limit = 14
        lhs = LHSHashTable(saved_hash, lhs_search_limit, [g_const_dict['zeta'](3)]) 

        poly_search_domain = Zeta3Domain1(
            [(2,2), (1,1), (1,17), (1,5)], # an coefs
            (-16,-1) # bn coef
            ) 

        enumerator = EfficentGCFEnumerator(
            lhs,
            poly_search_domain,
            [g_const_dict['zeta'](3)],
            lhs_search_limit
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

    def test_MITM_api4(self):
        saved_hash = 'pi_lhs_dept20'
        lhs_search_limit = 20
        lhs = LHSHashTable(
            saved_hash,
            lhs_search_limit,
            [g_const_dict['pi']]) 

        poly_search_domain = CartesianProductPolyDomain(
            1, [-13,13], # an deg 2, coefs ranging from -5 to 5
            2, [-11,11]) # bn deg 2, coefs ranging from -5 to 5

        # create an enumeator to iter thought the poly domain and compare it to the 
        # lhs table
        enumerator = EfficentGCFEnumerator(
            lhs,
            poly_search_domain,
            [g_const_dict['pi']],
            lhs_search_limit
            )

        results = enumerator.full_execution()
        self.assertEqual(len(results), 92)

    # def test_MITM_api5(self):
        
    #     cmd = 'python main.py MITM_RF -lhs_constant catalan -num_of_cores 2 -lhs_search_limit 20 -poly_a_order 3' +\
    #           ' -poly_a_coefficient_max 7 -poly_b_order 4 -poly_b_coefficient_max 2 --integer_factorization_bn'
    #     cmd = cmd.split(' ')[2:]
    #     parser = main.init_parser()
    #     args = parser.parse_args(cmd)
    #     results = main.enumerate_over_gcf_main(args)
    #     print(results)
    #     self.assertEqual(len(results), 1)
    #     self.assertIn('\\frac{2}{-1 + 2 Catalan\\left(\\right)} = 3 - \\frac{6}{13 - \\frac{64}{29 - \\frac{270}{51 - \\frac{768}{79 - \\frac{1750}{..}}}}}', results)

if __name__ == '__main__':
    unittest.main()
