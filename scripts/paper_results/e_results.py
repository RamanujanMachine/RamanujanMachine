from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator
from ramanujan.poly_domains.CartesianProductPolyDomain import CartesianProductPolyDomain
from ramanujan.constants import g_const_dict

"""
This script enumerates GCFs for e.

Both an and bn take the following form:
Xn = x0 * n^2 + x1 * n + x3
All coefs can range from -5 to 5.
All coefs are independent from one another and all combinations will be generated.
"""

# create a LHS table for e
saved_hash = 'e_lhs_dept5_db'
lhs_search_limit = 5
lhs = LHSHashTable(
    saved_hash,
    lhs_search_limit,
    [g_const_dict['e']]) 

# define the poly domain
poly_search_domain = CartesianProductPolyDomain(
    1, [-5, 5],
    1, [-5, 5])

# create an enumerator that iters thought the poly domain and compare GCFs to the lhs table
enumerator = EfficientGCFEnumerator(
    lhs,
    poly_search_domain,
    [g_const_dict['e']]
    )

results = enumerator.full_execution()
print("{} results found!".format(len(results)))
enumerator.print_results(results, 'unicode', convergence_rate=False)
