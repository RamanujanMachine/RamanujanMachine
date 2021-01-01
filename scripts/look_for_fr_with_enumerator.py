from FREnumerator import * 
from FixedBNPolyDomain import *

general_poly_domain = FixedBNPolyDomain(2, (-20, 20), (8,0,0,0,0)) 
# -(2n-1)n^3
# 1 + n + 3n^2
# 1 + 3n + 3n^2
'''
	poly_domain
	1. don't allow leading coef to be 0
	2. leading an >=1
	3. choose specific bn/an
'''

enumerator = FREnumerator(
	[], 
	1, 
	'', 
	general_poly_domain)


results = enumerator.find_initial_hits(general_poly_domain)
print (results)