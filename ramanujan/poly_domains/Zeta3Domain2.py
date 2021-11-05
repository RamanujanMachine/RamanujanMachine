from .CartesianProductPolyDomain import CartesianProductPolyDomain
import numpy as np


class Zeta3Domain2(CartesianProductPolyDomain):
	"""
	This domain iterates polynomials from this kind:

	a(n) = x0*n^3 + x0*(n+1)^3 + x1(2n+1)
	b(n) = -(x2**2)*n^6

	Note that all zeta3 results shown in the paper can be written using this scheme.

	It is suggested to keep x1,x2<0, but this is not enforced by this class
	"""

	def __init__(self, a_coefs_ranges=((0, 0),), b_coef_range=(0, 0), *args, **kwargs):
		# a_coef_range and b_coef_range are given blank values. they are initialized again afterwards
		super().__init__(a_deg=2, b_deg=1, a_coef_range=[0, 0], b_coef_range=[0, 0], *args, **kwargs)
		self.a_coef_range = a_coefs_ranges
		self.b_coef_range = [b_coef_range]

		self._setup_metadata()

	@staticmethod
	def get_calculation_method():
		def an_iterator(a_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield a_coefs[0] * ((i + 1) ** 3 + i ** 3) + a_coefs[1] * (2 * i + 1)

		def bn_iterator(b_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield -(b_coefs[0] ** 2) * (i ** 6)

		return an_iterator, bn_iterator

	def get_an_degree(self, an_coefs):
		return 3

	def get_bn_degree(self, bn_coefs):
		return 6

	def filter_gcfs(self, an_coefs, bn_coefs):
		# This scheme is very similar to other zeta schemes, so some duplications may occur
		a_leading_coef = an_coefs[0] * 2
		if -1 * (bn_coefs[0] ** 2) * 4 < -1 * (a_leading_coef ** 2):
			return False

		if self.use_strict_convergence_cond and -1 * (bn_coefs[0] ** 2) * 4 == -1 * (a_leading_coef ** 2):
			return False

		# discarding expansions
		if np.gcd.reduce(an_coefs + bn_coefs) != 1:
			return False

		return True
