import os
import sys
import argparse
import pickle
from time import time
from functools import partial
import sympy
from enumerate_over_gcf import multi_core_enumeration_wrapper
from enumerate_over_signed_rcf import esma_search_wrapper
import series_generators
import lhs_generators
import constants    # declares constants as sympy Singeltons, "not" used is intended

g_const_dict = {
    'zeta': sympy.zeta,
    'e': sympy.E,
    'pi': sympy.pi,
    'catalan': sympy.Catalan,
    'golden_ratio': sympy.GoldenRatio,
    'khinchin': sympy.S.Khinchin,
    'euler-mascheroni': sympy.EulerGamma,
    'pi-acosh_2': sympy.pi * sympy.acosh(2)
}


def get_custom_generator(generator_name, args):
    """
    Hackish custom generators. create a new one and add it here to include in API.
    :param generator_name: name of custom generator
    :param args: program args
    :return: generator function , number of free coefficients
    """
    if generator_name is None:
        return None, None
    elif generator_name == 'zeta_bn':
        return partial(series_generators.create_zeta_bn_series, args.function_value * 2), 2
    elif generator_name == 'zeta3_an':
        return series_generators.zeta3_an_generator, 2
    elif generator_name == 'zeta5_an':
        return series_generators.zeta5_an_generator, 3
    elif generator_name == 'catalan_bn':
        return series_generators.catalan_bn_generator, 2
    elif generator_name == 'polynomial_shift1':
        return series_generators.create_series_from_compact_poly_with_shift1, None
    elif generator_name == 'polynomial_shift2n1':
        return series_generators.create_series_from_compact_poly_with_shift2n1, None
    else:
        print("unknown custom series generator!")
        sys.exit(1)


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



def get_constant(const_name, args):
    """
    for function constants
    """
    if const_name == 'zeta':
        return g_const_dict[const_name](args.function_value)
    else:
        return g_const_dict[const_name]


def get_hash_filename(sympy_consts, args):
    name = ''
    for sym_const in sympy_consts:
        name += f'{str(sym_const)}_'
    name += f'{args.lhs_search_limit}_hash.p'
    name = name.replace('*', '_mul_')
    name = name.replace('/', '_div_')
    name = name.replace(' ', '')
    return name


# Initialize the argument parser that accepts inputs from the end user
def init_parser():
    parser = argparse.ArgumentParser(
        description='TODO',
        usage='''main.py <enumeration_type> [<args>]
        
Currently the optional enumeration types are:
    enumerate_over_gcf    enumerates over different RHS permutation by using 2 series generators for {an} and {bn}
    ESMA                  builds LHS enumeration, and conducts searches using the ESMA algorithm.
        ''')

    subparsers = parser.add_subparsers(help='enumeration type to run')

    gcf_parser = subparsers.add_parser('enumerate_over_gcf')
    gcf_parser.set_defaults(which='enumerate_over_gcf')
    gcf_parser.add_argument('-lhs_constant', choices=g_const_dict.keys(), nargs='+',
                            help='constants to search for - initializing the LHS hash table')
    gcf_parser.add_argument('-function_value', type=int,
                            help='Which value of the function are we assessing \
                                 (assuming LHS constant takes an arguments)')
    gcf_parser.add_argument('-lhs_search_limit', type=int,
                            help='The limit for the LHS coefficients')
    gcf_parser.add_argument('-num_of_cores', type=int,
                            help='The number of cores to run on', default=1)
    gcf_parser.add_argument('-poly_a_order', type=int,
                            help='the number of free coefficients for {a_n} series')
    gcf_parser.add_argument('-poly_a_coefficient_max', type=int,
                            help='The maximum value for the coefficients of the {a_n} polynomial')
    gcf_parser.add_argument('-poly_b_order', type=int,
                            help='the number of free coefficients for {b_n} series')
    gcf_parser.add_argument('-poly_b_coefficient_max', type=int,
                            help='The maximum value for the coefficients of the {b_n} polynomial')
    gcf_parser.add_argument('-custom_generator_an', type=str,
                            help='(optional) custom generator for {a_n} series. if defined, poly_a_order is ignored')
    gcf_parser.add_argument('-custom_generator_bn', type=str,
                            help='(optional) custom generator for {a_n} series. if defined, poly_b_order is ignored')

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


def enumerate_over_gcf_main(args):
    # Need to handle all the empty cases here, or have default values for them
    if not args.lhs_constant:
        print("You must input a parameter for the lhs_constant. Please run 'python main.py enumerate_over_gcf --help' to for more details about the required parameters.")
        return

    # same path to hash_tables no matter what
    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    # {an} series generator
    an_generator, poly_a_order = get_custom_generator(args.custom_generator_an, args)
    if (an_generator is None) or (poly_a_order is None):
        print('default series generator is chosen: a_n=n(n(...(c[0]*n + c[1]) + c[2]) + ...) + c[k]')
        poly_a_order = args.poly_a_order

    # {bn} series generator
    bn_generator, poly_b_order = get_custom_generator(args.custom_generator_bn, args)
    if (bn_generator is None) or (poly_b_order is None):
        print('default series generator is chosen: b_n=n(n(...(c[0]*n +ca[1]) + c[2]) + ...) + c[k]')
        poly_b_order = args.poly_b_order

    # constants for LHS
    sympy_consts = []
    for const in args.lhs_constant:
        sympy_consts.append(get_constant(const, args))
    hash_table_filename = get_hash_filename(sympy_consts, args)

    # Runs the enumeration wrapper
    final_results = multi_core_enumeration_wrapper(
        sym_constant=sympy_consts,  # constant to run on
        lhs_search_limit=args.lhs_search_limit,
        poly_a=[[i for i in range(args.poly_a_coefficient_max+1)]] * poly_a_order,  # a_n polynomial coefficients
        poly_b=[[i for i in range(args.poly_b_coefficient_max+1)]] * poly_b_order,  # b_n polynomial coefficients
        num_cores=args.num_of_cores,  # number of cores to run on
        manual_splits_size=None,  # use naive tiling
        saved_hash=os.path.join('hash_tables', hash_table_filename),  # if this doesn't exist, it will be created.
        create_an_series=an_generator,
        create_bn_series=bn_generator
    )

    with open('tmp_results', 'wb') as file:
        pickle.dump(final_results, file)
    return final_results


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
    elif args.which == 'enumerate_over_gcf':
        enumerate_over_gcf_main(args)
    elif args.which == 'ESMA':
        enumerate_over_signed_rcf_main(args)


if __name__ == '__main__':
    main()