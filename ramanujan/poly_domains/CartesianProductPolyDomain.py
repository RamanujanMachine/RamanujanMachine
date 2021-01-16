from .AbstractPolyDomains import * 
from ..utils.series_generators import iter_series_items_from_compact_poly
from itertools import product 

class CartesianProductPolyDomain(AbstractPolyDomains):
	"""
		This poly domain will generate all combinations for a(n) and b(n)
	"""
	def __init__(self, a_deg, a_coef_range, b_deg, b_coef_range, 
		*args, **kwargs):
		self.a_deg = a_deg
		# expanding the range to a different range for each coef
		# allows us to use the same functions for decedent classes
		self.a_coef_range = [a_coef_range for i in range(a_deg + 1)]
		self.b_deg = b_deg
		self.b_coef_range = [b_coef_range for i in range(b_deg + 1)]

		self.an_length = self.get_an_length()
		self.bn_length = self.get_bn_length()
		self.num_iterations = self.an_length * self.bn_length

		self.an_domain_range, self.bn_domain_range = self.dump_domain_ranges()

		super().__init__(*args, **kwargs)

	def _range_size(self, coef_range):
		return coef_range[1] - coef_range[0] + 1

	def domain_size_by_var_ranges(self, var_ranges):
		size = 1
		for var_range in var_ranges:
			size *= self._range_size(var_range)
		return  size

	# def expand_var_ranges_to_domain(self, var_ranges):
	# 	return [ \
	# 		[i for i in range(coef[0], coef[1] + 1)] \
	# 			for coef in var_ranges
	# 	]

	def get_an_length(self):
		return self.domain_size_by_var_ranges(self.a_coef_range)

	def get_bn_length(self):
		return self.domain_size_by_var_ranges(self.b_coef_range)

	def get_calculation_method(self):
		# both an and bn are regular compact polys
		return iter_series_items_from_compact_poly, \
			iter_series_items_from_compact_poly

	def dump_domain_ranges(self):
		an_domain = [ \
			[i for i in range(coef[0], coef[1] + 1)] \
				for coef in self.a_coef_range
		]

		bn_domain = [ \
			[i for i in range(coef[0], coef[1] + 1)] \
				for coef in self.b_coef_range
		]

		return an_domain, bn_domain

	def iter_polys(self, primary_looped_domain):
		an_domain, bn_domain = self.dump_domain_ranges()

		if primary_looped_domain == 'a':
			a_coef_iter = product(*an_domain)
			for a_coef in a_coef_iter:
				b_coef_iter = product(*bn_domain)
				for b_coef in b_coef_iter:
					yield a_coef, b_coef
		else:
			b_coef_iter = product(*bn_domain)
			for b_coef in b_coef_iter:
				a_coef_iter = product(*an_domain)
				for a_coef in a_coef_iter:
					yield a_coef, b_coef

	def get_a_coef_iterator(self):
		return product(*self.an_domain_range)
	
	def get_b_coef_iterator(self):
		return product(*self.bn_domain_range)
	
	def get_individual_polys_generators(self):
		# for backwards competability.
		an_domain, bn_domain = self.dump_domain_ranges()

		return product(*an_domain), product(*bn_domain)

