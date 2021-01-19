import os
import sys
import argparse
import pickle
from time import time
import sympy
from enumerate_over_signed_rcf import esma_search_wrapper
import ramanujan.utils.series_generators as series_generators
import lhs_generators
import ramanujan.constants  # declares constants as sympy Singeltons, "not" used is intended
from ramanujan.constants import g_const_dict


def get_lhs_generator(generator_name, args):
    """
    LHS enumeration generators. create a new one and add it here to include in API.
    If program is given an out_dir argument, will pickle and save the enumeration locally for later use.
    :param generator_name: name of custom generator
    :param args: program args
    :return: generator function , number of free coefficients
    """
    if generator_name == 'biased_monoms':
        if args.poly_deg <= 0 or len(args.coeff_lim) != 2 or args.coeff_lim[0] <= 0 or args.coeff_lim[1] <= 0:
            print("Degree and limits must be positive. Expecting two coefficient limits.")
            raise AttributeError
        return lhs_generators.create_biased_monoms(args.poly_deg, args.coeff_lim[0], args.coeff_lim[1])
    if generator_name == 'standard':
        if args.poly_deg <= 0 or not isinstance(args.coeff_lim, int) or args.coeff_lim <= 0:
            print("Degree and limits must be positive. Expecting a single coefficient limit.")
            raise AttributeError
        return lhs_generators.create_standard_lhs(args.poly_deg, args.coeff_lim, args.out_dir, do_print=(not args.no_print))


# Initialize the argument parser that accepts inputs from the end user
def init_parser():
    parser = argparse.ArgumentParser(
        description='TODO',
        usage='''main.py <enumeration_type> [<args>]
        
Currently the you can use only this enumerator:
    ESMA                  builds LHS enumeration, and conducts searches using the ESMA algorithm.
        ''')

    """ ESMA PARSER """

    subparsers = parser.add_subparsers(help='ESMA')

    srcf_parser = subparsers.add_parser('ESMA')

    srcf_parser.set_defaults(which='ESMA')
    srcf_parser.add_argument('-out_dir', type=str,
                             help='Directory to output binary results or file path for existing generic enumeration')
    srcf_parser.add_argument('-mode', choices=['build', 'search'],
                             help='Build LHS enumerations or find conjectures')
    # Search-only arguments:
    srcf_parser.add_argument('-constant', choices=g_const_dict.keys(), nargs='?', default=None, const=None,
                             help='constant to search for')
    srcf_parser.add_argument('-cycle_range', type=int, nargs='+', default=None, const=None,
                             help='Shortest and longest period-length for the signed sequences to search')
    srcf_parser.add_argument('-depth', type=int, nargs='?', default=None, const=None,
                             help='In case depth needs to be changed (if insufficient precision error repeats)')
    srcf_parser.add_argument('-no_print', action='store_true')

    # Dual-purpose arguments:
    srcf_parser.add_argument('-lhs', type=str, nargs='?', default=None, const=None,
                             help='Name of LHS pattern to build, or path of LHS enumeration to use')
    srcf_parser.add_argument('-poly_deg', type=int, nargs='?', default=None, const=None,
                             help='Maximum degree of numerator and denominator polynomials in the LHS function')
    srcf_parser.add_argument('-min_deg', type=int, nargs='?', default=None, const=None,
                             help='Minimum degree from which to search. Useful for breaking up searches.')
    srcf_parser.add_argument('-coeff_lim', type=int, nargs='?', default=None, const=None,
                             help='Maximum absolute value for the coefficients of the LHS function')
    return parser

def enumerate_over_signed_rcf_main(args):
    if args.mode == 'build':
        if args.lhs is None:
            print("When building LHS enumeration, a structure must be specified")
            raise ValueError
        print('Building an LHS enumeration:')
        if args.out_dir is not None:
            if os.path.exists(args.out_dir):
                print("File under given name already exists. Choose different name for output file.")
                return
            else:
                print('Saving the pickled enumeration to ' + str(args.out_dir))
        lhs = get_lhs_generator(args.lhs, args)
        return lhs
    if args.mode == 'search':
        print('Running a search for conjectures using ESMA algorithm:')
        if args.lhs is not None:
            with open(args.lhs, 'rb') as f:
                print("Starting to load existing LHS enumeration:")
                strt = time()
                custom_lhs = pickle.load(f)
                print("Loaded {} LHS variations. Took {} sec".format(len(custom_lhs), time() - strt))
        else:
            custom_lhs = None
            if min(args.poly_deg, args.coeff_lim) <= 0:
                print("poly_deg and coeff_lim must be positive integers.")
                raise ValueError
        if len(args.cycle_range) != 2 or args.cycle_range[0] > args.cycle_range[1]:
            print("Cycle range must be two non-decreasing  positive integers.")
            raise ValueError
        results, _ = esma_search_wrapper(constant=g_const_dict[args.constant],
                                         custom_enum=custom_lhs,
                                         poly_deg=args.poly_deg,
                                         coeff_lim=args.coeff_lim,
                                         cycle_range=args.cycle_range,
                                         min_deg=args.min_deg,
                                         depth=args.depth,
                                         out_dir=args.out_dir,
                                         do_print=(not args.no_print))
        return results


def main():
    # Initializes the argument parser to receive inputs from the user
    parser = init_parser()

    args = parser.parse_args()

    if len(sys.argv) == 1:
        print("You must input the running parameters to run the Ramanujan Machine. Please run 'python main.py \
         enumerate_over_gcf  --help' to for more details about the required parameters.")
    elif args.which == 'ESMA':
        enumerate_over_signed_rcf_main(args)


if __name__ == '__main__':
    main()
