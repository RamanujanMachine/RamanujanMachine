import sympy
from ramanujan.poly_domains.Zeta5Domain import *
from ramanujan.HugeLHSHashTable import *
from ramanujan.LHSHashTable import *
from ramanujan.enumerators.GuessFromFREnumerator import *
from ramanujan.enumerators.HugeDomainRelativeGCFEnumerator import *
from ramanujan.enumerators.EfficentGCFEnumerator import *

saved_hash = 'e.lhs.dept5.db'
const_e = [lambdify((), sympy.E, modules="mpmath")()]
lhs_search_limit = 5
poly_search_domain = CartesianProductPolyDomain(
	2, [-5,5], # an deg 2, coefs ranging from -5 to 5
	2, [-5,5]) # bn deg 2, coefs ranging from -5 to 5


lhs = LHSHashTable(
    saved_hash,
    lhs_search_limit,
    [sympy.E], # constant
    10**-10) 


with mpmath.workdps(10000):
    lhs = LHSHashTable(
        saved_hash,
        lhs_search_limit,
        const_e,  # constant
        10**-10)  # length of key

enumerator = EfficentGCFEnumerator(
	lhs,
	poly_search_domain,
	const_e,
	lhs_search_limit
	)

enumerator = EfficentGCFEnumerator(
	lhs,
	poly_search_domain,
	[sympy.E],
	lhs_search_limit
	)

results = enumerator.find_initial_hits(poly_search_domain)
better_results = enumerator.refine_results(results)

with mpmath.workdps(10000):

python main.py 
	MITM_RF 
	-lhs_constant e V
	-num_of_cores 1 V
	-lhs_search_limit 5 V
	-poly_a_order 2 V
	-poly_a_coefficient_max 5 V
	-poly_b_order 2  V
	-poly_b_coefficient_max 5 V
