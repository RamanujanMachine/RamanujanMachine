import unittest
import main
import os
import pickle
import sympy
from lhs_generators import create_standard_lhs


class APITests(unittest.TestCase):

    def test_ESMA_api1(self): # Test full enumeration and search configuration including saving binaries. Might ake a little longer
        cmd = 'ESMA -out_dir ./tmp -mode search -constant e -cycle_range 2 2 -depth 105 -poly_deg 1' + \
              ' -coeff_lim 2 -no_print'
        cmd = cmd.split(' ')
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        results = main.enumerate_over_signed_rcf_main(args)
        self.assertEqual(len(results), 13)
        adjusted = [[res[0], res[1], list(res[3])] for res in results]
        self.assertIn([(sympy.E / (sympy.E - 1)), [1, -1], [1, 0, -2, 0, 1]], adjusted)
        self.assertIn([-1 + sympy.E, [-1, 1], [1, 0, -2, 0, 1]], adjusted)
        print('Search results are as expected.')
        files_there = os.path.exists('./tmp/res_list_0') and os.path.exists('./tmp/recurring_by_value_0')
        self.assertTrue(files_there)
        os.remove('./tmp/res_list_0')
        os.remove('./tmp/recurring_by_value_0')
        os.rmdir('./tmp')
        print("Successfully removed result output files.")

    def test_ESMA_api2(self): # Test standard build configuration.
        cmd = 'ESMA -out_dir ./tmp -mode build -lhs standard -poly_deg 1 -coeff_lim 1 -no_print'
        cmd = cmd.split(' ')
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        lhs = main.enumerate_over_signed_rcf_main(args)
        print('Creating enumeration not through API to compare:')
        self.assertEqual(lhs, create_standard_lhs(poly_deg=1, coefficients_limit=1, do_print=(not args.no_print)))
        print("Identical enumerations.")
        file_there = os.path.exists('./tmp')
        self.assertTrue(file_there)
        os.remove('./tmp')
        print('Successfuly removed output file')

    def test_ESMA_api3(self): # Test search using an existing enumeration configuration.
        print('Creating and saving a temporary generic LHS enumeration.')
        custom_enum = create_standard_lhs(poly_deg=1, coefficients_limit=2, do_print=False)
        path = './tmp'
        with open(path, 'wb') as file:
            pickle.dump(custom_enum, file)
        print('Calling using API:')
        cmd = 'ESMA, -mode, search, -constant, e, -cycle_range, 2, 2, -lhs, ./tmp, -no_print'
        cmd = cmd.split(', ')
        parser = main.init_parser()
        args = parser.parse_args(cmd)
        print('Searching using generic LHS')
        results = main.enumerate_over_signed_rcf_main(args)
        os.remove(path)
        print('Deleted temporary generic LHS enumeration from disk')
        self.assertEqual(len(results), 13)
        adjusted = [[res[0], res[1], list(res[3])] for res in results]
        self.assertIn([(sympy.E / (sympy.E - 1)), [1, -1], [1, 0, -2, 0, 1]], adjusted)
        self.assertIn([(sympy.E / (-2 + sympy.E)), [1, 1], [1, 0, 0, -1, 0, 0, -1, 0, 0, 1]], adjusted)
        print('Search results are as expected.')


if __name__ == '__main__':
    unittest.main()
