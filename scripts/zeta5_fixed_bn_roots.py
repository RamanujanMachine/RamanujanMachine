from ramanujan.enumerators.FREnumerator import FREnumerator
from ramanujan.poly_domains.Zeta5DomainFromSpecificBnRoots import Zeta5DomainFromSpecificBnRoots
from ramanujan.constants import g_const_dict
from ramanujan.multiprocess_enumeration import multiprocess_enumeration

"""
"""


def main():
    # define the poly domain

    poly_search_domain = Zeta5DomainFromSpecificBnRoots(
        [(1, 1), (-0, 0), (-50, 50), (-50, 50)],
        bn_roots=[(2, 1), (2, -1),
                  (2, 3), (2, -3),
                  (2, 5), (2, -5),
                  (2, 7), (2, -7),
                  (2, 9), (2, -9)
            ])

    # enum = FREnumerator(poly_search_domain,
    #     [g_const_dict['zeta'](5), g_const_dict['zeta'](3)])

    # results = enum.full_execution()

    results = multiprocess_enumeration(
        FREnumerator,
        None,  # No LHS in FREnumerator
        poly_search_domain,
        [g_const_dict['zeta'](5), g_const_dict['zeta'](3)],
        4)

    print("{} results found!".format(len(results)))
    for r in results:
        print(r)


if __name__ == '__main__':
    main()
