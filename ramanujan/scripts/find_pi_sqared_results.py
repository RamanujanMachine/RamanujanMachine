import sympy

from MonomExponentialGCFEnumerator import MonomExponentialGCFEnumerator
from CartesianProductPolyDomain import CartesianProductPolyDomain
from RelativeGCFEnumerator import RelativeGCFEnumerator
from PiSqaredDomain1 import PiSqaredDomain1
from PiSqaredDomain2 import PiSqaredDomain2

sym_constant = [sympy.pi ** 2]
lhs_search_limit = 30
saved_hash = 'hash_tables/pi_mul__mul_2_30_hash.p'
poly_domains_generator = CartesianProductPolyDomain(2, [1,3], 4, [-5,5])


pi_sqared_dom2 = PiSqaredDomain2([
	(3,3), (-20,20), (-20,20),  # a coefs
 	(1,6), (-10, 10), (1,6), (-10, 10)]) # b coefs

enumerator = RelativeGCFEnumerator(sym_constant, lhs_search_limit, saved_hash, pi_sqared_dom2)

res = enumerator._first_enumeration(pi_sqared_dom2, True)
print('finished_first_enum')
res = enumerator._refine_results(res, True)




# pis_dom = PiSqaredDomain1([(1,20), (-20,20), (-20,20), (-20, 20)])

# enumerator = RelativeGCFEnumerator(sym_constant, lhs_search_limit, saved_hash, pis_dom)

# res = enumerator._first_enumeration(pis_dom, True)
# print('finished_first_enum')
# res = enumerator._refine_results(res, True)