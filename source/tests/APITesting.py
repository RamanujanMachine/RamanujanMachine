import unittest
import main
import os


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


if __name__ == '__main__':
    unittest.main()
