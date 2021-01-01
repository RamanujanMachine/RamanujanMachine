from .CartesianProductPolyDomain import * 

class PiSqaredDomain1(CartesianProductPolyDomain):
	'''
	This domain iters polynomials from this kind:
	a(n) = x0*n^2 + x0*n + x1
	b(n) = x2*n^4 + x3*n^3

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
			a_deg = 2, # deg reffers to degree of freedom, not poly_deg
			b_deg = 2,
			a_coef_range = None, # no use, overridden
			b_coef_range = None
			,*args, **kwargs)

	def get_calculation_method(self):
		def an_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*i*(i+1) + free_vars[1]

		def bn_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				yield free_vars[0]*(i**4) + free_vars[1]*(i**3)

		return an_iterator, bn_iterator

	def get_an_length(self):
		return self._range_size(self.coef_ranges[0]) * \
			 self._range_size(self.coef_ranges[1])

	def get_bn_length(self):
		return self._range_size(self.coef_ranges[2]) * \
			 self._range_size(self.coef_ranges[3])

	def dump_domain_ranges(self):
		an_domain = [ \
			[i for i in range(self.coef_ranges[0][0], self.coef_ranges[0][1])],
			[i for i in range(self.coef_ranges[1][0], self.coef_ranges[1][1])],
		]
		bn_domain = [ \
			[i for i in range(self.coef_ranges[2][0], self.coef_ranges[2][1])],
			[i for i in range(self.coef_ranges[3][0], self.coef_ranges[3][1])],
		]

		return an_domain, bn_domain