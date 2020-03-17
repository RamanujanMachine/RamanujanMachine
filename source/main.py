import os
import sys
import argparse
import pickle
from functools import partial
import sympy
from enumerate_over_gcf import multi_core_enumeration_wrapper
import enumerate_over_signed_rcf
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
    'euler-mascheroni': sympy.gamma,
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
        if args.poly_deg <= 0 or len(args.coeff_lim) != 1 or args.coeff_lim <= 0:
            print("Degree and limits must be positive. Expecting a single coefficient limit.")
            raise AttributeError
        return lhs_generators.create_std_lhs(args.poly_deg, args.coeff_limit, args.out_dir)



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

        ''')

    subparsers = parser.add_subparsers(help='enumeration type to run')

    gcf_parser = subparsers.add_parser('mitm_rf')
    gcf_parser.set_defaults(which='mitm_rf')
    gcf_parser.add_argument('-LHS_constant', choices=g_const_dict.keys(), nargs='+',
                            help='constants to search for - initializing the LHS hash table')
    gcf_parser.add_argument('--function_value', type=int,
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
    gcf_parser.add_argument('--custom_generator_an', type=str,
                            help='(optional) custom generator for {a_n} series. if defined, poly_a_order is ignored')
    gcf_parser.add_argument('--custom_generator_bn', type=str,
                            help='(optional) custom generator for {a_n} series. if defined, poly_b_order is ignored')

    srcf_parser = subparsers.add_parser('massey')
    srcf_parser.set_defaults(which='massey_rama')
    srcf_parser.add_argument('-out_dir', type=str,
                             help='Directory to output binary results or generic enumeration')
    srcf_parser.add_argument('-mode', choices=['build', 'search'],
                             help='Build LHS enumerations or find conjectures')
    #Search-only arguments:
    srcf_parser.add_argument('-constant', choices=g_const_dict.keys(), nargs=1, default=None, const=None,
                             help='constant to search for')
    srcf_parser.add_argument('-cycle_range', type=int, nargs='?', default=None, const=None,
                             help='Shortest and longest period-length for the signed sequences to search')
    srcf_parser.add_argument('-depth', type=int, nargs='?', default=None, const=None,
                             help='In case depth needs to be changed (if insufficient precision error repeats)')

    #Dual-purpose arguments:
    srcf_parser.add_argument('-lhs', type=str,
                             help='Name of LHS pattern to build, or path of LHS enumeration to use')
    srcf_parser.add_argument('-poly_deg', type=int, nargs='?', default=None, const=None,
                                 help='Maximum degree of numerator and denominator polynomials in the LHS function')
    srcf_parser.add_argument('-min_deg', type=int, nargs='?', default=None, const=None,
                             help='Minimum degree from which to search. Useful for breaking up searches.')
    srcf_parser.add_argument('-coeff_lim', type=int, nargs='?', default=None, const=None,
                             help='Maximum absolute value for the coefficients of the LHS function')
    return parser


def enumerate_over_gcf_main(args):
    # same path to hash_tables no matter what
    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    # {an} series generator
    an_generator, poly_a_order = get_custom_generator(args.custom_generator_an, args)
    if (an_generator is None) or (poly_a_order is None):
        poly_a_order = args.poly_a_order

    # {bn} series generator
    bn_generator, poly_b_order = get_custom_generator(args.custom_generator_bn, args)
    if (bn_generator is None) or (poly_b_order is None):
        poly_b_order = args.poly_b_order

    # constants for LHS
    sympy_consts = []
    for const in args.LHS_constant:
        sympy_consts.append(get_constant(const, args))
    hash_table_filename = get_hash_filename(sympy_consts, args)

    # Runs the enumeration wrapper
    final_results = multi_core_enumeration_wrapper(
        sym_constant=sympy_consts,  # constant to run on
        lhs_search_limit=args.lhs_search_limit,
        poly_a=[[i for i in range(args.poly_a_coefficient_max)]] * poly_a_order,  # a_n polynomial coefficients
        poly_b=[[i for i in range(args.poly_b_coefficient_max)]] * poly_b_order,  # b_n polynomial coefficients
        num_cores=args.num_of_cores,  # number of cores to run on
        manual_splits_size=None,  # use naive tiling
        saved_hash=os.path.join('hash_tables', hash_table_filename),  # if existing
        create_an_series=an_generator,  # use default
        create_bn_series=bn_generator
    )

    with open('tmp_results', 'wb') as file:
        pickle.dump(final_results, file)
    return final_results


def enumerate_over_signed_rcf_main(args):
    if args.mode == 'build':
        if args.lhs is None:
            print("When building LHS enumeration, a structure must be specified")
            raise AttributeError
        lhs = get_lhs_generator(args.LHS_structure, args)
    if args.mode == 'search':
        if len(args.cycle_range) != 2 or min(args.poly_deg, args.coeff_lim, min(args.cycle_range)) <= 0 or \
                args.cycle_range[0] >= args.cycle_range[1]:
            print("Degree, limits and range must be positive. Cycle range must be two increasing  positive integers.")
            raise AttributeError
        res_list, dup_dict = enumerate_over_signed_rcf.search_wrapper(args.constant, args.lhs, args.poly_deg,
                                                                      args.coeff_lim, args.cycle_range,
                                                                      args.min_deg, args.depth)
        if args.out_dir:
            path = args.out_dir
        else:
            path = 'tmp'
        if not os.path.exists(path):
            os.makedirs(path)
        res = '/'.join([path, 'dups_by_value'])
        dup = '/'.join([path, 'res_list'])
        with open(res, 'wb') as f:
            pickle.dump(res_list, f)
        with open(dup, 'wb') as f:
            pickle.dump(dup_dict, f)



def main():
    # Initializes the argument parser to receive inputs from the user
    parser = init_parser()
    if len(sys.argv) == 1:  # run from editor
        args = parser.parse_args(['enumerate_over_gcf'])
    else:
        args = parser.parse_args()
    if args.which == 'mitm_rf':
        enumerate_over_gcf_main(args)
    elif args.which == 'massey_rama':
        enumerate_over_signed_rcf_main(args)


if __name__ == '__main__':
    main()
