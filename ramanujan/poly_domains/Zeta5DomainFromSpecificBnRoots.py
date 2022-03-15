from .CartesianProductPolyDomain import CartesianProductPolyDomain
import numpy as np


class Zeta5DomainFromSpecificBnRoots(CartesianProductPolyDomain):
	"""
	an is a degenerated 5th deg, of this form:
		an = n^5 + (n+1)^5 + a0*n^3 + a1*n^2 + a3*n + a4

	bn is defined by constant roots. Makes things easier.
	"""

	def __init__(self, a_coefs_ranges=((0, 0),), bn_roots=[(0, 0)], *args, **kwargs):
		"""
		Bn roots is of this form - [(root1_num, root1_denom), (root2_num, root2_denom), ...]
		Bn is of 10th degree. If not enough roots are given, the rest will be 0 (root = (1,0)). 
		"""
		# a_coef_range and b_coef_range are given blank values. they are initialized again afterwards
		super().__init__(a_deg=4, b_deg=0, a_coef_range=[0, 0], b_coef_range=[0, 0], *args, **kwargs)
		self.a_coef_range = a_coefs_ranges
		self.b_coef_range = []
		if len(bn_roots) != 10:
			bn_roots += [(1,0)] * (10 - len(bn_roots))
		self.bn_roots = bn_roots

		# To keep desc = 1, we need to choose the same leading coef for an, and sqrt(lead_coef(bn))
		# Since the roots are constant, this calculation may happen only once here.
		self.an_leading_coef = 1
		for root in bn_roots:
			self.an_leading_coef *= root[0]
		self.an_leading_coef = np.sqrt(self.an_leading_coef)

		# make sure that bn roots create an integer coef for an
		if not self.an_leading_coef.is_integer():
			raise Exception('Non-integer coefficient for an. Roots choise is not valid.')
		self.an_leading_coef = int(self.an_leading_coef)

		self._setup_metadata()

	def get_calculation_method(self):
		def an_iterator(a_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield self.an_leading_coef * (
						i**5 + (i+1)**5 + \
						a_coefs[0]*i**3 + a_coefs[1]*i**2 + a_coefs[2]*i +a_coefs[3]
					)

		def bn_iterator(b_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				next_item = 1 
				for root in self.bn_roots:
					next_item *= (root[0]*i + root[1])
				yield -next_item

		return an_iterator, bn_iterator

	def get_an_degree(self, an_coefs):
		return 5

	def get_bn_degree(self, bn_coefs):
		return 10

	def filter_gcfs(self, an_coefs, bn_coefs):
		# maybe add test if bn roots are also roots for an
		return True
