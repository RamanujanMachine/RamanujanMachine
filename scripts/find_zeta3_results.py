import sympy

from MonomExponentialGCFEnumerator import MonomExponentialGCFEnumerator
from CartesianProductPolyDomain import CartesianProductPolyDomain
from RelativeGCFEnumeratorDebug import RelativeGCFEnumeratorDebug
from Zeta3Domain1 import Zeta3Domain1

sym_constant = [sympy.zeta(3)]
lhs_search_limit = 10
saved_hash = 'hash_tables/zeta(3)_10_has1h.p'

zeta_dom = Zeta3Domain1(
	[(2,10), (1,10), (-200,200), (-200, 200), # an 
	(-30,-2)] #bn
	)
enumerator = RelativeGCFEnumeratorDebug(sym_constant, lhs_search_limit, saved_hash, zeta_dom)

res = enumerator.find_initial_hits(zeta_dom, True)
print('finished_first_enum')
refined_res = enumerator.refine_results(res)
refined_res

# pis_dom = PiSqaredDomain1([(1,20), (-20,20), (-20,20), (-20, 20)])

# enumerator = RelativeGCFEnumerator(sym_constant, lhs_search_limit, saved_hash, pis_dom)

# res = enumerator._first_enumeration(pis_dom, True)
# print('finished_first_enum')
# res = enumerator._refine_results(res, True)

from HugeLHSHashTable import *
weird_zeta3 = lambdify((), sympy.zeta(3), modules="mpmath")()
with mpmath.workdps(50):
	HugeLHSHashTable(
		'zeta_dept_20', 
		100,
		[weird_zeta3], 
		30,
		skip_bloom_load=True)


from HugeLHSHashTable import *
weird_zeta3 = lambdify((), sympy.zeta(3), modules="mpmath")()
with mpmath.workdps(50):
	lhs = HugeLHSHashTable(
		'zeta_dept_20', 
		20,
		[weird_zeta3], 
		10)


import sympy

from CartesianProductPolyDomain import CartesianProductPolyDomain
from HugeDomainRelativeGCFEnumerator import HugeDomainRelativeGCFEnumerator
from Zeta3Domain1 import Zeta3Domain1

sym_constant = [sympy.zeta(3)]
lhs_search_limit = 10
saved_hash = 'zeta_dept100_30_partitions'

zeta_dom = Zeta3Domain1(
	[(2,2), (1,1), (1,10), (10, 20), # an 
	(-4,-1)] #bn
	)
enumerator = HugeDomainRelativeGCFEnumerator(
	sym_constant, 
	lhs_search_limit, 
	saved_hash, 
	zeta_dom)

res = enumerator.find_initial_hits(zeta_dom, True)
print('finished_first_enum')
refined_res = enumerator.refine_results(res)
refined_res


with mpmath.workdps(50):
    constants = [const() for const in self.constants_generator]
    lhs = HugeLHSHashTable(
        saved_hash, 
        self.lhs_limit,
        constants, 
        g_N_initial_key_length)

import sympy

from CartesianProductPolyDomain import CartesianProductPolyDomain
from HugeDomainRelativeGCFEnumerator import HugeDomainRelativeGCFEnumerator
from Zeta3Domain1 import Zeta3Domain1

sym_constant = [sympy.zeta(3)]
lhs_search_limit = 10
saved_hash = 'zeta_dept100_30_partitions'

zeta_dom = Zeta3Domain1(
	[(2,2), (1,1), (1,300), (-300, 300), # an 
	(-50,-1)] #bn
	)
enumerator = HugeDomainRelativeGCFEnumerator(
	sym_constant, 
	lhs_search_limit, 
	saved_hash, 
	zeta_dom)

res = enumerator.find_initial_hits(zeta_dom, True)
print('finished_first_enum')
refined_res = enumerator.refine_results(res)
refined_res


from Zeta3Domain1 import *
from HugeLHSHashTable import *
weird_zeta3 = lambdify((), sympy.zeta(3), modules="mpmath")()
with mpmath.workdps(50):
	lhs = HugeLHSHashTable(
		'zeta_dept100_30_partitions', 
		10,
		[weird_zeta3], 
		10,
		skip_bloom_load=True)



zeta_dom = Zeta3Domain1(
	[(-100,100), (-300,300), (-300,300), # an 
	(-1,-50)] #bn
	)