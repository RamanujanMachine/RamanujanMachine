from ramanujan.enumerators.FREnumerator import FREnumerator
from ramanujan.poly_domains.Zeta3Domain2 import Zeta3Domain2
from ramanujan.constants import g_const_dict
from ramanujan.multiprocess_enumeration import multiprocess_enumeration


"""
This script enumerates GCFs for zeta(3) using the FR algorithm.
It uses a specific family of series that tends to generate conjectures for zeta(3)
The family is defined under ramanujan.poly_domains.Zeta3Domain3

It uses the following series:
a(n) = x0(n**3 + (n+1)**3) + x1(2n+1)
b(n) = x2*n^6
"""


def main():
    # define the poly domain
    poly_search_domain = Zeta3Domain2(
        [(1, 3), (-20, 20)],
        (1, 2))

    # create an enumerator to iter thought the poly domain and compare it to the lhs table
    enumerator = FREnumerator(
        poly_search_domain,
        [g_const_dict['zeta'](3)]
    )

    # results = enumerator.full_execution()
    results = multiprocess_enumeration(
        FREnumerator,
        None, # No LHS in FREnumerator
        poly_search_domain,
        [g_const_dict['zeta'](3)],
        4)

    print("{} results found!".format(len(results)))
    for r in results:
        print(r)

if __name__ == '__main__':
    main()