from .CartesianProductPolyDomain import CartesianProductPolyDomain 


class Zeta3DomainWithRatC(CartesianProductPolyDomain):
	"""
	A sub-set of Zeta3Domain2 is an infinite family of GCFs with FR, that converges to permutation of zeta(3):
		an = n^3 (n+1)^3 + 2c(c+1)(2n+1)
		bn = -n^6
	When c can be any rational number. However, we want to use only integer coefficients. We can achieve this by
	inflation:
		c <- c + x/y
		an = n^3 + (n+1)^3 + 2(c+x/y)(c+x/y+1)(2n+1)
	inflate by y^2, and get:
		an = y^2(n^3 + (n+1)^3) + 2(cy+x)(cy+x+y)(2n+1)
		bn = -y^4 * b^6
	"""
	def __init__(self, a_coefs_ranges=((0, 0),), b_coef_range=(0, 0), *args, **kwargs):
		# a_coef_range and b_coef_range are given blank values. they are initialized again afterwards
		super().__init__(a_deg=3, b_deg=1, a_coef_range=[0, 0], b_coef_range=[0, 0], *args, **kwargs)
		self.a_coef_range = a_coefs_ranges
		self.b_coef_range = [b_coef_range]

		self._setup_metadata()

	@staticmethod
	def get_calculation_method():
		def an_iterator(a_coefs, max_runs, start_n=1):
			# to match notation in class doc - a_coefs = [c, x, y]
			for i in range(start_n, max_runs):
				yield (a_coefs[2]**2) * ((i+1)**3 + i**3) + \
					2 * (a_coefs[0]*a_coefs[2] + a_coefs[1]) * \
					(a_coefs[0]*a_coefs[2] + a_coefs[1] + a_coefs[2]) * (2*i+1)

		def bn_iterator(b_coefs, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield -(b_coefs[0]**4)*(i**6)

		return an_iterator, bn_iterator

	def get_an_degree(self, an_coefs):
		return 3

	def get_bn_degree(self, bn_coefs):
		return 6

	def filter_gcfs(self, an_coefs, bn_coefs):
		"""
		This class is not truly a cartesian domain - a and b coefficients share a value (y).
		This means that GCFs with different values for a's y and b's y will be generated. We filter them out here.
		"""
		# Require same y for a and b
		if an_coefs[2] != bn_coefs[0]:
			return False
		# We don't use x/y >= 1, this is a duplicate of using another c
		if an_coefs[1] >= an_coefs[2]:
			return False

		return True
