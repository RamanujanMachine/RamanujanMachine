import os
import sys
import argparse
import pickle
from functools import partial
import sympy
from enumerate_over_gcf import multi_core_enumeration_wrapper
import series_generators
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
        poly_a_order = args.poly_a_order

    # {bn} series generator
    bn_generator, poly_b_order = get_custom_generator(args.custom_generator_bn, args)
    if (bn_generator is None) or (poly_b_order is None):
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


def main():
    # Initializes the argument parser to receive inputs from the user
    parser = init_parser()

    if len(sys.argv) == 1:  # run from editor
        args = parser.parse_args(['enumerate_over_gcf'])
    else:
        args = parser.parse_args()

    if len(sys.argv) == 1:
        print("You must input the running parameters to run the Ramanujan Machine. Please run 'python main.py enumerate_over_gcf --help' to for more details about the required parameters.")
    elif args.which == 'enumerate_over_gcf':
        enumerate_over_gcf_main(args)


if __name__ == '__main__':
    main()
