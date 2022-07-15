from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.constants import g_const_dict

"""
This script enumerates GCFs for zeta(3).
It uses a specific family of series that tends to generate conjectures for zeta(3)
The family is defined under ramanujan.poly_domains.Zeta3Domain1

It uses the following series:
a(n) = (x0*n + x1)(x2*n*(n + 1) + x3)
b(n) = x4*n^6
We're setting x0 = 2 and x1 = 1.
x2 and x3 are set to range from -50 to 50
x4 is set to range from -10 to -1

Also, Zeta3Domain1 will check for a convergence condition for the GCF, as described in http://www.ramanujanmachine.com/
"""

# create a LHS table for zeta
saved_hash = "zeta3.lhs.dept20.db"
lhs_search_limit = 20
lhs = LHSHashTable(saved_hash, lhs_search_limit, [g_const_dict["zeta"](3)])

# define the poly domain
poly_search_domain = Zeta3Domain1([(2, 2), (1, 1), (-50, 50), (-50, 50)], (-10, -1))

# create an enumerator to iter thought the poly domain and compare it to the lhs table
enumerator = EfficientGCFEnumerator(lhs, poly_search_domain, [g_const_dict["zeta"](3)])

results = enumerator.full_execution()
print("{} results found!".format(len(results)))
enumerator.print_results(results, "unicode", convergence_rate=False)
