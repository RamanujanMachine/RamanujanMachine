from ramanujan.poly_domains.Zeta3Domain1 import *
from ramanujan.HugeLHSHashTable import *
from ramanujan.enumerators.GuessFromFREnumerator import *

weird_zeta3 = lambdify((), sympy.zeta(3), modules="mpmath")()
with mpmath.workdps(50):
	lhs = HugeLHSHashTable(
		'zeta_dept100_30_partitions', 
		1,
		[weird_zeta3], 
		10,
		skip_bloom_load=True)

# zeta_dom = Zeta3Domain1(
# 	[(-100,100), (-300,300), (-300,300), # an 
# 	(-1,-1)] #bn
# 	)

zeta_dom = Zeta3Domain1(
	[(2,2), (-10,10), (-10,10), (-10,10), # an 
	(-1,-1)] #bn
	)

enumerator = GuessFromFREnumerator(
	[], 
	1, 
	'', 
	zeta_dom)

results = enumerator.find_initial_hits(zeta_dom)
print (results)
better_results = enumerator.refine_results(results)
print (better_results)