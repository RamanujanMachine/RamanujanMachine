from AbstractPolyDomains import * 
from series_generators import iter_series_items_from_compact_poly
from itertools import product 

class FROrientedPolyDomain(AbstractPolyDomains):
	"""
		Will iterate over bns  by monomes, and skip ans that are an expansion of a previus bn

		If we'll generate monoms in the simplest way possible, we'll generate through plenty
		of identical bns - n(n+1)(n+1) = (n+1)(n+1)n = (n+1)n(n+1)

		hence, the number of UNIQUE combinations is [(number of possible monomes) choose i] for i 
		ranging from 1 to b_deg	
	"""
	def __init__(self, a_deg, a_coef_range, b_deg, b_n_coef_range, b_free_coef_range, 
		*args, **kwargs):
		self.a_deg = a_deg 
		self.a_coef_range = a_coef_range
		self.b_deg = b_deg
		self.b_n_coef_range = b_n_coef_range
		self.b_free_coef_range = b_free_coef_range

		self.an_length = self.get_an_length()
		self.bn_length = self.get_bn_length()
		self.num_iterations = self.an_length * self.bn_length

		super().__init__(*args, **kwargs)

	def get_an_length(self):
		return (self.a_deg + 1) ** (self.a_coef_range[1] - self.a_coef_range[0] + 1)

	def get_bn_length(self):
		number_of_monoms = (self.b_n_coef_range[1] - self.b_n_coef_range[0]) * \
			(self.b_free_coef_range[1] - self.b_free_coef_range[0])
		return self.b_deg ** number_of_monoms

	def get_calculation_method(self):
		def bn_iterator(free_vars, max_runs, start_n=1):
			for i in range(start_n, max_runs):
				result = 1
				for j in b_deg:	
					result *= (free_vars[j]*i + free_vars[j+1])
				yield result
		
		return iter_series_items_from_compact_poly, bn_iterator

	@staticmethod
	def iter_monoms(n_coef_range, free_coef_range):
		for n_coef in range(n_coef_range[0], n_coef_range[1]):
			if n_coef == 0:
				continue
			for free_coef in range(free_coef_range[0], free_coef_range[1]):
				yield (n_coef, free_coef)

	@staticmethod
	def recursive_poly_generator(monoms_left, poly, n_coef_range, free_coef_range):
		all_monoms = [i for i in FROrientedPolyDomain.iter_monoms(n_coef_range, free_coef_range)]
		for monom in all_monoms:
			current_poly = poly + [monom]
			if monoms_left == 1:
				sorted_poly = list(current_poly)
				sorted_poly.sort() 

				# only one premutation should be sorted, this way we discard duplicates
				if sorted_poly == current_poly:
					unpacked_poly = []
					for m in current_poly:
						unpacked_poly.append(m[0])
						unpacked_poly.append(m[1])
					yield unpacked_poly
			else:
				yield from recursive_poly_generator(monoms_left - 1, current_poly, n_coef_range, free_coef_range)

	def check_for_convegence(self, an_coefs, bn_coefs):
		b_leading_coef = 1
		for i in range(0, bn_coefs, 2):
			b_leading_coef *= bn_coefs[i]

		# checking for >= as well as >, might be overkill
		return b_leading_coef * 4 >= -1 * (an_coefs[0]**2)

	def iter_polys(self, primary_looped_domain):
		an_domain = [[i for i in range(self.a_coef_range[0], self.a_coef_range[1]+1)] \
			for i in range(self.a_deg + 1)]

		if primary_looped_domain == 'a':
			a_coef_iter = product(*an_domain)
			for a_coef in a_coef_iter:
				b_coef_iter = recursive_poly_generator(self.b_deg, [], self.b_n_coef_range, self.b_free_coef_range)
				for b_coef in b_coef_iter:
					if self.check_for_convege(a_coef, b_coef):
						yield a_coef, b_coef
		else:
			b_coef_iter = recursive_poly_generator(self.b_deg, [], self.b_n_coef_range, self.b_free_coef_range)
			for b_coef in b_coef_iter:
				a_coef_iter = product(*an_domain)
				for a_coef in a_coef_iter:
					if self.check_for_convege(a_coef, b_coef):
						yield a_coef, b_coef

	def dump_domain_ranges(self):
		pass