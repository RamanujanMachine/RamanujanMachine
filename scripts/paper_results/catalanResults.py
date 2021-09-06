from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator
from ramanujan.poly_domains.CatalanDomain import CatalanDomain
from ramanujan.constants import g_const_dict

"""

"""

# create a LHS table for catalan
saved_hash = 'catalan_lhs_dept14.db'
lhs_search_limit = 8
lhs = LHSHashTable(
    saved_hash,
    lhs_search_limit,
    [g_const_dict['catalan'], g_const_dict['pi-acosh_2']])

# define the poly domain
poly_search_domain = CatalanDomain(
    a_coefs_ranges=(-14, 14),
    poly_an_degree=2,
    b_coefs_ranges=((-5,5), (-5,5))
    )

# create an enumerator to iter thought the poly domain and compare it to the lhs table
enumerator = EfficientGCFEnumerator(
    lhs,
    poly_search_domain,
    [g_const_dict['catalan'], g_const_dict['pi-acosh_2']]
)

results = enumerator.full_execution()
print("{} results found!".format(len(results)))
enumerator.print_results(results, 'unicode', convergence_rate=False)
