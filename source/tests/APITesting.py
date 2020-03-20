import unittest
import main
import os
from sympy import E as e
from enumerate_over_signed_rcf import SignedRcfEnumeration
from lhs_generators import create_std_lhs
import pickle

class APITests(unittest.TestCase):
    def test_api1(self):
        cmd = ['enumerate_over_gcf', '-lhs_constant', 'e', '-num_of_cores', '1', '-lhs_search_limit', ' 5',
               '-poly_a_order', ' 2', '-poly_a_coefficient_max', '4', '-poly_b_order', ' 2', '-poly_b_coefficient_max', '4']
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        results = main.enumerate_over_gcf_main(args)
        print(results)
        self.assertEqual(len(results), 17)
        self.assertIn('\\frac{1 + e}{-1 + e} = 2 + \\frac{1}{6 + \\frac{1}{10 + \\frac{1}{14 + \\frac{1}{18 + \\frac{1}{..}}}}}',
                      results)
        self.assertIn('\\frac{1}{-2 + e} = 1 + \\frac{1}{2 + \\frac{2}{3 + \\frac{3}{4 + \\frac{4}{5 + \\frac{5}{..}}}}}',
                      results)

    def test_api2(self):
        cmd = ['enumerate_over_gcf', '-lhs_constant', 'zeta', '-function_value', '3', '-num_of_cores', '2',
               '-lhs_search_limit', '14', '-poly_a_order', '3', '-poly_a_coefficient_max', '19',
               '-poly_b_order', '3', '-poly_b_coefficient_max', '19',
               '-custom_generator_an', 'zeta3_an', '-custom_generator_bn', 'zeta_bn']
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        results = main.enumerate_over_gcf_main(args)
        print(results)
        self.assertEqual(len(results), 3)
        self.assertIn('\\frac{8}{7 \\zeta\\left(3\\right)} = 1 - \\frac{1}{21 - \\frac{64}{95 - \\frac{729}{259 - \\frac{4096}{549 - \\frac{15625}{..}}}}}',
                      results)
        self.assertIn('\\frac{12}{7 \\zeta\\left(3\\right)} = 2 - \\frac{16}{36 - \\frac{1024}{160 - \\frac{11664}{434 - \\frac{65536}{918 - \\frac{250000}{..}}}}}',
                      results)
        self.assertIn('\\frac{6}{\\zeta\\left(3\\right)} = 5 - \\frac{1}{117 - \\frac{64}{535 - \\frac{729}{1463 - \\frac{4096}{3105 - \\frac{15625}{..}}}}}',
                      results)

    def test_api3(self):    # this one take a few minutes
        cmd = ['enumerate_over_gcf', '-lhs_constant', 'catalan', 'pi-acosh_2', '-num_of_cores', '1',
               '-lhs_search_limit', '8', '-poly_a_order', '3', '-poly_a_coefficient_max', '14',
               '-poly_b_order', '1', '-poly_b_coefficient_max', '5',
               '-custom_generator_bn', 'catalan_bn']
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        results = main.enumerate_over_gcf_main(args)
        print(results)
        self.assertEqual(len(results), 1)
        self.assertIn('\\frac{6}{- \\pi \\operatorname{acosh}{\\left(2 \\right)} + 8 Catalan\\left(\\right)} = 2 - \\frac{2}{19 - \\frac{108}{56 - \\frac{750}{113 - \\frac{2744}{190 - \\frac{7290}{..}}}}}',
                      results)

    def test_api4(self):  # this one take a few minutes
        cmd = ['enumerate_over_gcf', '-lhs_constant', 'pi', '-num_of_cores', '2',
               '-lhs_search_limit', '20', '-poly_a_order', '2', '-poly_a_coefficient_max', '13',
               '-poly_b_order', '3', '-poly_b_coefficient_max', '11',
               '-custom_generator_bn', 'polynomial_shift1']
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        results = main.enumerate_over_gcf_main(args)
        print(results)
        self.assertEqual(len(results), 20)

    def test_api5(self): # Test full enumeration and search configuration.
        cmd = ['ESMA', '-mode', 'search', '-constant', 'e',
               '-cycle_range', '2', '2', '-depth', '105', '-poly_deg', '1', '-coeff_lim', '2', '-no_print']
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        results = main.enumerate_over_signed_rcf_main(args)
        self.assertEqual(len(results), 13)
        adjusted = [[res[0], res[1], list(res[3])] for res in results]
        self.assertIn([(e / (e - 1)), [1, -1], [1, 0, -2, 0, 1]], adjusted)
        self.assertIn([-1 + e, [-1, 1], [1, 0, -2, 0, 1]], adjusted)

    def test_api6(self): # Test standard build configuration.
        cmd = ['ESMA', '-out_dir', './tmp', '-mode', 'build',
               '-lhs', 'standard', '-poly_deg', '1', '-coeff_lim', '1', '-no_print']
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        lhs = main.enumerate_over_signed_rcf_main(args)
        print('Creating enumeration not through API to compare:')
        self.assertEqual(lhs, create_std_lhs(poly_deg=1, coefficients_limit=1, no_print=args.no_print))
        print("Identical enumerations.")
        file_there = os.path.exists('./tmp')
        self.assertTrue(file_there)
        os.remove('./tmp')
        print('Successfuly removed output file')

    def test_api7(self): # Test search using an existing enumeration configuration.
        print('Creating and saving a temporary generic LHS enumeration.')
        custom_enum = create_std_lhs(poly_deg=1, coefficients_limit=2, no_print=True)
        path = './tmp'
        with open(path, 'wb') as file:
            pickle.dump(custom_enum, file)
        print('Calling using API:')
        cmd = ['ESMA', '-mode', 'search', '-constant', 'e',
               '-cycle_range', '2', '2', '-lhs', './tmp', '-no_print']
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        print('Searching using generic LHS')
        results = main.enumerate_over_signed_rcf_main(args)
        os.remove(path)
        print('Deleted temporary generic LHS enumeration from disk')
        self.assertEqual(len(results), 13)
        adjusted = [[res[0], res[1], list(res[3])] for res in results]
        self.assertIn([(e / (e - 1)), [1, -1], [1, 0, -2, 0, 1]], adjusted)
        self.assertIn([(e / (-2 + e)), [1, 1], [1, 0, 0, -1, 0, 0, -1, 0, 0, 1]], adjusted)
        print('Search results are as expected.')
if __name__ == '__main__':
    unittest.main()
