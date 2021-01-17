import sympy
from ramanujan.LHSHashTable import *
from ramanujan.enumerators.EfficentGCFEnumerator import *
from ramanujan.poly_domains.CartesianProductPolyDomain import *

# create a LHS table for e
saved_hash = 'e_lhs_dept5_db'
lhs_search_limit = 5
lhs = LHSHashTable(
    saved_hash,
    lhs_search_limit,
    [g_const_dict['e']]) 

# define a poly domain to iter through
poly_search_domain = CartesianProductPolyDomain(
    2, [-5,5], # an deg 2, coefs ranging from -5 to 5
    2, [-5,5]) # bn deg 2, coefs ranging from -5 to 5

# create an enumeator to iter thought the poly domain and compare it to the 
# lhs table
enumerator = EfficentGCFEnumerator(
    lhs,
    poly_search_domain,
    [g_const_dict['e']],
    lhs_search_limit
    )

results = enumerator.full_execution()