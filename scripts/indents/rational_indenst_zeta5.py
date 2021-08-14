from ramanujan.poly_domains.RationalIndents import *
from ramanujan.poly_domains.Zeta5Indents import *
from ramanujan.enumerators.IndentsEnumerator import IndentsEnumerator
from ramanujan.constants import g_const_dict
from ramanujan.multiprocess_enumeration import multiprocess_enumeration


def an(n, coeffs):
	return coeffs[0]*((n+1)**5 + n**5) + coeffs[1]*((n+1)**3 + n**3) + coeffs[2]*(2*n+1)


def bn(n,coeffs):
	return -(coeffs[0]**2)*(n**10)


def main():
	dom = Zeta5Indents(
		[
			[[1, 0, 0], [1]],
			[[1, 6, 0], [1]],
			[[1, 6, -4], [1]],
			[[1, 16, -4], [1]]
		],
		[an, bn],
		20,
		5, 10
		)

	results = multiprocess_enumeration(
		IndentsEnumerator,
		None,
		dom,
		[g_const_dict['zeta'](5)],
		4
		)

	print(results)

if __name__ == '__main__':
	main()