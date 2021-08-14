from ramanujan.poly_domains.RationalIndents import *
from ramanujan.enumerators.IndentsEnumerator import IndentsEnumerator
from ramanujan.constants import g_const_dict
from ramanujan.multiprocess_enumeration import multiprocess_enumeration
import json
import time

def an(n, coeffs):
	return coeffs[0]*((n+1)**3 + n**3) + coeffs[1]*(2*n+1)


def bn(n,coeffs):
	return -(coeffs[0]**2)*(n**6)


def main():
	dom = RationalIndents(
		[
			[[1, 4], [1]],
			[[1, 0], [1]],
			[[5, -3], [4]],
			[[17, -12], [1]],
			[[3, -2], [1]],
			[[1, 4], [1]],
			[[1, 12], [1]],
			[[1, 24], [1]],
			[[1, 15], [1]],
			[[1, 3], [1]]
		],
		[an, bn],
		20,
		3, 6
		)

	# en = IndentsEnumerator(dom, [g_const_dict['zeta'](3)])
	# results = en.full_execution()
	results = multiprocess_enumeration(
		IndentsEnumerator,
		None,
		dom,
		[g_const_dict['zeta'](3)],
		4
		)

	with open('zeta3_indents_' + time.strftime("%y%m%d_%H%M%S") + '.json', 'w') as f:
		json.dump([res._asdict() for res in results], f, indent=2)

if __name__ == '__main__':
	main()