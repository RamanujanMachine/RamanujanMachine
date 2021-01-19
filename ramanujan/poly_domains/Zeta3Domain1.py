from .CartesianProductPolyDomain import CartesianProductPolyDomain
from itertools import product


class Zeta3Domain1(CartesianProductPolyDomain):
	"""
	This domain iters polynomials from this kind:
	a(n) = (x0*n + x1)(x2*n*(n + 1) + x3)
	b(n) = x4*n^6

	where x0, x1, x2, x3, x4 are 5 freedom degrees

	to reduce search space and be similar to other results-
	keep x0 and x1 low
	keep x4 negative

	this is a decedent of CartesianProductPolyDomain since an and bn has no	particular relation

	NOTICE - Since every coeffiecnt is given explicitly, we do not enforce that the leading coef of an will always be
	positive. (See CartesianProductPolyDomain documnetation regarding an_leading_coef_positive for more information)
	"""
	def __init__(self, a_coefs_ranges, b_coef_range, *args, **kwargs):
		"""
		:param a_coefs_ranges: the range allowed for each coef from x0,x1,x2,x3
		in this format-
			[(x0_min, x0_max), ... ]
		:param b_coefs_ranges: b has only one coef, so this will hold (x4_min, x4_max)
		"""
		# a_coef_range and b_coef_range are given blank values. they are initialized again afterwards
		super().__init__(a_deg=4, b_deg=1, a_coef_range=[0, 0], b_coef_range=[0, 0], *args, **kwargs)
		self.a_coef_range = a_coefs_ranges
		self.b_coef_range = [b_coef_range]

		self.an_length = self.get_an_length()
		self.bn_length = self.get_bn_length()
		self.num_iterations = self.an_length * self.bn_length

		self.an_domain_range, self.bn_domain_range = self.dump_domain_ranges()

	def get_calculation_method(self):
		def an_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield (free_vars[0]*i + free_vars[1])*(free_vars[2]*i*(i+1) + free_vars[3])

		def bn_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*(i**6)

		return an_iterator, bn_iterator

	@staticmethod
	def get_poly_an_degree(an_coefs):
		deg = 3
		if an_coefs[0] == 0:
			deg -= 1
		if an_coefs[2] == 0:
			deg -= 2
		return deg

	@staticmethod
	def get_poly_bn_degree(bn_coefs):
		# bn_coefs is not used since the degree is always 0. Still accepting this variable for consistency
		return 6

	@staticmethod
	def get_poly_an_lead_coef(an_coefs):
		return an_coefs[0] * an_coefs[2]

	@staticmethod
	def get_poly_bn_lead_coef(bn_coefs):
		return bn_coefs[0]
	
	@staticmethod
	def check_for_convergence(an_coefs, bn_coefs):
		# see Ramanujan paper for convergence condition on balanced an & bn degrees
		a_leading_coef = an_coefs[0] * an_coefs[2]

		# checking for >= as well as >, might be overkill
		return bn_coefs[0] * 4 >= -1 * (a_leading_coef**2)

	def iter_polys(self, primary_looped_domain):
		an_domain, bn_domain = self.dump_domain_ranges()

		if primary_looped_domain == 'a':
			a_coef_iter = product(*an_domain)
			for a_coef in a_coef_iter:
				b_coef_iter = product(*bn_domain)
				for b_coef in b_coef_iter:
					if self.check_for_convergence(a_coef, b_coef):
						yield a_coef, b_coef
		else:
			b_coef_iter = product(*bn_domain)
			for b_coef in b_coef_iter:
				a_coef_iter = product(*an_domain)
				for a_coef in a_coef_iter:
					if self.check_for_convergence(a_coef, b_coef):
						yield a_coef, b_coef
