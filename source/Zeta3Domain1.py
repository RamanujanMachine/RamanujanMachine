from CartesianProductPolyDomain import * 

class Zeta3Domain1(CartesianProductPolyDomain):
	'''
	This domain iters polynomials from this kind:
	a(n) = (x0*n + x1)(x2*n*(n + 1) + x3)
	b(n) = x4*n^6

	where x0, x1, x2, x3, x4 are 5 freedom degrees

	to reduce search space and be similar to other results- 
	keep x0 and x1 low
	keep x4 negative

	this is a decendent of CartesianProductPolyDomain since an and bn has no
	praticular relation
	'''

	def __init__(self, coef_ranges, *args, **kwargs):
		'''
		coef_ranges - the range allowd for each coef from x0,x1,x2,x3
		in this format-
			[(x0_min, x0_max), ... ]
		'''
		self.coef_ranges = coef_ranges

		super().__init__(
			a_deg = 4, # deg reffers to degree of freedom, not poly_deg
			b_deg = 1,
			a_coef_range = None, # no use, overridden
			b_coef_range = None
			,*args, **kwargs)

	def get_calculation_method(self):
		def an_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield (free_vars[0]*i + free_vars[1])*(free_vars[2]*i*(i+1) + free_vars[3])

		def bn_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*(i**6)

		return an_iterator, bn_iterator

	def get_an_length(self):
		return self.domain_size_by_var_ranges(self.coef_ranges[:4])

	def get_bn_length(self):
		return self.domain_size_by_var_ranges(self.coef_ranges[4:])

	def dump_domain_ranges(self):
		return self.expand_var_ranges_to_domain(self.coef_ranges[:4]), \
			self.expand_var_ranges_to_domain(self.coef_ranges[4:])

	# allows user to get actual data regarding the polynomals he got
	@classmethod
	def get_poly_an_degree(an_coefs):
		deg = 3
		if an_coefs[0] == 0:
			deg -= 1
		if an_coefs[2] ==0:
			deg -= 2
		return deg

	@classmethod
	def get_poly_bn_degree(bn_coefs):
		return 6

	@classmethod
	def get_poly_an_lead_coef(an_coefs):
		return an_coefs[0] * an_coefs[2]

	@classmethod
	def get_poly_bn_lead_coef(bn_coefs):
		return bn_coefs[0]

	
	def get_individual_polys_generators(self):
		# skips throught non conveging examples
		an_domain, bn_domain = self.dump_domain_ranges()

		return product(*an_domain), product(*bn_domain)

