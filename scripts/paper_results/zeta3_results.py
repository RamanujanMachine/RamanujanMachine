import sympy
from ramanujan.LHSHashTable import *
from ramanujan.enumerators.EfficentGCFEnumerator import *
from ramanujan.poly_domains.Zeta3Domain1 import *

# create a LHS table for zeta
saved_hash = 'zeta3.lhs.dept20.db'
lhs_search_limit = 20
lhs = LHSHashTable(
    saved_hash,
    lhs_search_limit,
    [g_const_dict['zeta'](3)]) 

# define a poly domain to iter through
poly_search_domain = Zeta3Domain1(
	[(2,2), (1,1), (-50,50), (-50,50)], # an coefs
	(-10,-1) # bn coef
	) 

# create an enumeator to iter thought the poly domain and compare it to the 
# lhs table
enumerator = EfficentGCFEnumerator(
	lhs,
	poly_search_domain,
	[g_const_dict['zeta'](3)],
	lhs_search_limit
	)

results = enumerator.full_execution()