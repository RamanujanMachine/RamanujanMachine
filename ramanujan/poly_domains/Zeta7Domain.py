from .CartesianProductPolyDomain import CartesianProductPolyDomain
import numpy as np


class Zeta7Domain(CartesianProductPolyDomain):
	"""
	Following the same idea proposed in Zeta5Domain, We use this structure as a domain for Zeta7:
		an = x0(n^7 + (n+1)^7) + x1(n^5 + (n+1)^5) + x2(n^3 + (n+1)^3) + x3(2n+1)
		bn = -x4^2 * n^14
	"""

	def __init__(self, a_coefs_ranges=((0,0),), b_coef_range=(0,0), *args, **kwargs):
		# a_coef_range and b_coef_range are given blank values. they are initialized again afterwards
		super().__init__(a_deg=4, b_deg=1, a_coef_range=[0, 0], b_coef_range=[0, 0], *args, **kwargs)
		self.a_coef_range = a_coefs_ranges
		self.b_coef_range = [b_coef_range]

		self._setup_metadata()

	def get_calculation_method(self):
		def an_iterator(a_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield a_coefs[0]*((i+1)**7 + i**7) + a_coefs[1]*((i+1)**5 + i**5) + \
					a_coefs[2]*(i**3 + (i+1)**3) + a_coefs[3]*(2*i+1)

		def bn_iterator(b_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield -(b_coefs[0]**2)*(i**14)

		return an_iterator, bn_iterator

	def get_an_degree(self, an_coefs):
		return 7

	def get_bn_degree(self, bn_coefs):
		return 14

	def filter_gcfs(self, an_coefs, bn_coefs):
		# This scheme is very simular to other zeta schemes, so some duplications may occure
		a_leading_coef = an_coefs[0] * 2
		if -1 * (bn_coefs[0]**2) * 4 < -1 * (a_leading_coef**2):
			return False

		if self.use_strict_convergence_cond and -1 * (bn_coefs[0]**2) * 4 == -1 * (a_leading_coef**2):
			return False

		# discarding expansions
		if np.gcd.reduce(an_coefs + bn_coefs) != 1:
			return False

		return True
