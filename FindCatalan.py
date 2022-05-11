# Imports
from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.constants import g_const_dict

from ramanujan.poly_domains.CartesianProductPolyDomain import CartesianProductPolyDomain

from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator

# Setup search
saved_hash = 'G_lhs_dept10'
lhs_coefs_range = 10
lhs = LHSHashTable(saved_hash,
lhs_coefs_range, [g_const_dict['catalan']])

poly_search_domain = CartesianProductPolyDomain(
4, [-5, 5], # an coefs
3, [-5, 5]) # bn coefs

enumerator = EfficientGCFEnumerator(
    lhs,
    poly_search_domain,
    [g_const_dict['catalan']]
    )

# Start execution
results = enumerator.full_execution()

# Display results
print("------------------------------------ R E S U L T S ------------------------------------")
enumerator.print_results(results)


