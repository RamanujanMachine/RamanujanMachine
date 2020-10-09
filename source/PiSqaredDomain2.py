from CartesianProductPolyDomain import * 

class PiSqaredDomain2(CartesianProductPolyDomain):
	'''
	This domain iters polynomials from this kind:
	a(n) = x0*n^2 + x1*n + x2
	b(n) = -n*(n+2)*(x3*n + x4)*(x5*n + x6)

	where x0, x1, x2, x3, x4, x5, x6 are 7 freedom degrees

	to reduce search space- 
	keep x0 at 3 or 5 ish
	keep x3, x5 low

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
			a_deg = 2, # deg reffers to degree of freedom, not poly_deg
			b_deg = 4,
			a_coef_range = None, # no use, overridden
			b_coef_range = None
			,*args, **kwargs)

	def get_calculation_method(self):
		def an_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*(i**2) + free_vars[1]*i + free_vars[2]

		def bn_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield -1*i*(i+2)*(free_vars[0]*i + free_vars[1])*(free_vars[2]*i + free_vars[3])

		return an_iterator, bn_iterator

	def get_an_length(self):
		return self.domain_size_by_var_ranges(self.coef_ranges[:3])

	def get_bn_length(self):
		return self.domain_size_by_var_ranges(self.coef_ranges[3:])

	def dump_domain_ranges(self):
		return self.expand_var_ranges_to_domain(self.coef_ranges[:3]), \
			self.expand_var_ranges_to_domain(self.coef_ranges[3:])