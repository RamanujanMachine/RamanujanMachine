from .CartesianProductPolyDomain import * 

class Zeta5Domain(CartesianProductPolyDomain):
	'''
	It apprears that zeta3 domain is something like:
		n^3 + (n+1)^3 + an + b
	so, it looks like n^2 is missing. Lets try to apply the same logic to zeta5!

	a(n) = x0*n^5 + x0*(n+1)^5 + x1*n^3 + x2*n^2 + x3*n + x4
	b(n) = x5*n^10

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
			a_deg = 5, # deg reffers to degree of freedom, not poly_deg
			b_deg = 1,
			a_coef_range = None, # no use, overridden
			b_coef_range = None
			,*args, **kwargs)

	def get_calculation_method(self):
		def an_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*( (i+1)**5 + i**5 ) + \
					free_vars[1] * (i**3) + free_vars[2] * (i**2) + \
					free_vars[3] * i + free_vars[4]

		def bn_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*(i**10)

		return an_iterator, bn_iterator

	def get_an_length(self):
		return self.domain_size_by_var_ranges(self.coef_ranges[:5])

	def get_bn_length(self):
		return self.domain_size_by_var_ranges(self.coef_ranges[5:])

	def dump_domain_ranges(self):
		return self.expand_var_ranges_to_domain(self.coef_ranges[:5]), \
			self.expand_var_ranges_to_domain(self.coef_ranges[5:])
	
	def get_individual_polys_generators(self):
		# skips throught non conveging examples
		an_domain, bn_domain = self.dump_domain_ranges()

		return product(*an_domain), product(*bn_domain)

