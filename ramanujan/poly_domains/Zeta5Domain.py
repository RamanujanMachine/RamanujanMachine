from .CartesianProductPolyDomain import * 
import numpy as np

class Zeta5Domain(CartesianProductPolyDomain):
	'''
	It apprears that zeta3 domain is something like:
		a(n) = x0(n^3 + (n+1)^3) + x1(2n+1)
		b(n) = x2*n^3
	so, it looks like n^2 is missing. Lets try to apply the same logic to zeta5!
		x0(n^5 + (n+1)^5) + x1(n^3 + (n+1)^3) + x2(2n+1)
		b(n) = -(x3**2)*n^10
	'''

	def __init__(self, a_coefs_ranges=((0,0),), b_coef_range=(0,0), *args, **kwargs):
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

	@staticmethod
	def get_an_degree(an_coefs):
		return 5

	@staticmethod
	def get_bn_degree(bn_coefs):
		return 10

	@staticmethod
	def get_poly_an_lead_coef(an_coefs):
		return an_coefs[0]*2

	@staticmethod
	def get_poly_bn_lead_coef(bn_coefs):
		return bn_coefs[0]

	def filter_gcfs(self, an_coefs, bn_coefs):
		# see Ramanujan paper for convergence condition on balanced an & bn degrees
		a_leading_coef = an_coefs[0] * 2
		if -1 * (bn_coefs[0]**2) * 4 < -1 * (a_leading_coef**2):
			return False

		if self.use_strict_convergence_cond and -1 * (bn_coefs[0]**2) * 4 == -1 * (a_leading_coef**2):
			return False

		# discarding expansions
		if np.gcd.reduce(an_coefs + bn_coefs) != 1:
			return False

		return True
