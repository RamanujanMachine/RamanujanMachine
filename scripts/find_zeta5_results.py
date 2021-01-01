import sympy
from ramanujan.poly_domains.Zeta5Domain import *
from ramanujan.HugeLHSHashTable import *
from ramanujan.enumerators.GuessFromFREnumerator import *
from ramanujan.enumerators.HugeDomainRelativeGCFEnumerator import *
from ramanujan.enumerators.FRAndCalcEnumerator import *

zeta5 = lambdify((), sympy.zeta(5), modules="mpmath")()


# option 1 - an 'skips' n^4
zeta_dom = Zeta5Domain(
	[(1,10), (-30,30), (-30,30), (-30,30), (-30,30), # an 
	(-1,-1)] #bn
	)

# option 2 - an is n^5 + (n+1)^5 + liniar factor
zeta_dom = Zeta5Domain(
	[(1,20), (0,0), (0,0), (-100,100), (-100,100), # an 
	(-1,-1)] #bn
	)


# try 1- 
# there is no LHS, just look for things that have FR under 
# this family
enumerator = GuessFromFREnumerator(
	[], 
	1, 
	'', 
	zeta_dom
	)

results = enumerator.find_initial_hits(zeta_dom)
print(results)


# try 2 - use clasic relative enumerator
with mpmath.workdps(20):
	lhs = HugeLHSHashTable(
 		'zeta5_dept50_20_partitions', 
 		50, # coefs range
 		[zeta5], 
 		20 # partitions
 		)


saved_hash = 'zeta5_dept50_20_partitions'

enumerator = HugeDomainRelativeGCFEnumerator(
	[zeta5], 
	50, 
	saved_hash, 
	zeta_dom)

res = enumerator.find_initial_hits(zeta_dom, True)
print('finished_first_enum')


# try 3 - use FR for first enumerator and HT for refinning
# try 1- 
# there is no LHS, just look for things that have FR under 
# this family
enumerator = FRAndCalcEnumerator(
	[zeta5], 
	50, 
	saved_hash, 
	zeta_dom)
