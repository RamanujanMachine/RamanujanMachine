from ramanujan.enumerators.FREnumerator import FREnumerator
from ramanujan.poly_domains.Zeta3DomainWithRatC import Zeta3DomainWithRatC
from ramanujan.constants import g_const_dict
from ramanujan.multiprocess_enumeration import multiprocess_enumeration


"""
This scripts demonstrate an infinite family of GCFs that converges to a permutation of zeta(3).
The family consists of those series:
    an = n^3 (n+1)^3 + 2c(c+1)(2n+1)
    bn = -n^6
While c can be any rational number. To describe c as a rational number we use
    c <- c + x / y
The GCF value will be:
    -2/polygamma(2,c+x/y). 
During the calculation we inflate (multiply by constant) the GCF, and the result might be this value multiplied by a 
factor.
For integer values (x=0) polygamma is a mobius transform of zeta(3). For rational number, it is a more complex 
expression and using PSLQ with zeta(3) will not recognize it, so this script doesn't give a definite value for them.
"""


def main():
    poly_search_domain = Zeta3DomainWithRatC(
        [(1, 5), (0, 5), (1, 5)],  # [c, x, y]
        (1, 5))
    results = multiprocess_enumeration(
        FREnumerator,
        None,  # No LHS in FREnumerator
        poly_search_domain,
        [g_const_dict['zeta'](3)],
        4)

    print("{} results found!".format(len(results)))
    for r in results:
        print(r)


if __name__ == '__main__':
    main()
