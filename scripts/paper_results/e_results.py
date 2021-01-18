from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator
from ramanujan.poly_domains.CartesianProductPolyDomain import CartesianProductPolyDomain
from ramanujan.constants import g_const_dict

# create a LHS table for e
saved_hash = 'e_lhs_dept5_db'
lhs_search_limit = 5
lhs = LHSHashTable(
    saved_hash,
    lhs_search_limit,
    [g_const_dict['e']]) 

# both an and bn take the following form:
# Xn = x0 * n^2 + x1 * n + x3
# all coefs can range from -5 to 5.
# all coefs are independent and all combinations will be generated.
poly_search_domain = CartesianProductPolyDomain(
    2, [-5, 5],
    2, [-5, 5])

# create an enumerator to iter thought the poly domain and compare it to the
# lhs table
enumerator = EfficientGCFEnumerator(
    lhs,
    poly_search_domain,
    [g_const_dict['e']],
    lhs_search_limit
    )

results = enumerator.full_execution()
