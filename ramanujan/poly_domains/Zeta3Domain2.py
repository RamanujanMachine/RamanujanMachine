from .CartesianProductPolyDomain import * 

class Zeta3Domain2(CartesianProductPolyDomain):
	'''
	This domain iters polynomials from this kind:
	a(n) = x0*n^3 + x0*(n+1)^3 + x1*x2(n+1) + x2
	b(n) = x3*n^6

	where x0, x1, x2, x3 are 4 freedom degrees

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
			a_deg = 3, # deg reffers to degree of freedom, not poly_deg
			b_deg = 1,
			a_coef_range = None, # no use, overridden
			b_coef_range = None
			,*args, **kwargs)

	def get_calculation_method(self):
		def an_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*( (i+1)**3 + i**3 ) + free_vars[1]*free_vars[2]*(i+1) + free_vars[2]

		def bn_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*(i**6)

		return an_iterator, bn_iterator

	def get_an_length(self):
		return self.domain_size_by_var_ranges(self.coef_ranges[:3])

	def get_bn_length(self):
		return self.domain_size_by_var_ranges(self.coef_ranges[3:])

	def dump_domain_ranges(self):
		return self.expand_var_ranges_to_domain(self.coef_ranges[:3]), \
			self.expand_var_ranges_to_domain(self.coef_ranges[3:])
	
	def get_individual_polys_generators(self):
		# skips throught non conveging examples
		an_domain, bn_domain = self.dump_domain_ranges()

		return product(*an_domain), product(*bn_domain)

