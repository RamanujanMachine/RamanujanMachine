from .CartesianProductPolyDomain import CartesianProductPolyDomain
import numpy as np


class Zeta5Domain(CartesianProductPolyDomain):
	"""
	This domain generalizes the idea build in zeta3 domain's to zeta5
	The zeta3 domain is:
		a(n) = x0(n^3 + (n+1)^3) + x1(2n+1)
		b(n) = x2*n^3
	It looks like n^2 is missing. In this domain we apply the same logic to zeta5!
		a(n) = x0(n^5 + (n+1)^5) + x1(n^3 + (n+1)^3) + x2(2n+1)
		b(n) = -(x3**2)*n^10
	We skip the even degrees.
	"""

	def __init__(self, a_coefs_ranges=((0, 0),), b_coef_range=(0, 0), *args, **kwargs):
		# a_coef_range and b_coef_range are given blank values. they are initialized again afterwards
		super().__init__(a_deg=3, b_deg=1, a_coef_range=[0, 0], b_coef_range=[0, 0], *args, **kwargs)
		self.a_coef_range = a_coefs_ranges
		self.b_coef_range = [b_coef_range]

		self._setup_metadata()

	def get_calculation_method(self):
		def an_iterator(a_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield a_coefs[0] * ((i+1)**5 + i**5) + a_coefs[1] * (i**3 + (i+1)**3) + a_coefs[2] * (2*i+1)

		def bn_iterator(b_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield -(b_coefs[0]**2)*(i**10)

		return an_iterator, bn_iterator

	def get_an_degree(self, an_coefs):
		return 5

	def get_bn_degree(self, bn_coefs):
		return 10

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
