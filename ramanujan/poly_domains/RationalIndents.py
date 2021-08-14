from .AbstractPolyDomains import AbstractPolyDomains
from ..utils.utils import iter_series_items_from_compact_poly
from sympy import var, Poly
from numpy import gcd, array_split
from copy import deepcopy


class RationalIndents(AbstractPolyDomains):
	def __init__(self, seeds, series_functions, max_denom, an_deg, bn_deg, starting_index=None, num_iters=None):
		# TODO - add MP support 
		self.an_deg = an_deg
		self.bn_deg = bn_deg
		self.seeds = seeds
		self.series_functions = series_functions
		"""
		def foo_an(n, coef_vec):
			return n*coef_vec[0] + n**2 * coef_vec[1] ... 
		"""
		self.max_denom = max_denom

		self.starting_index = starting_index if starting_index else 1

		if num_iters:
			self.num_iters = num_iters
		else:
			# This is a sum of the series 1, 2, 3, ...., max_denom - 1
			self.num_iters = max_denom * (max_denom - 1) / 2

		self.an_length = self.bn_length = self.num_iters

		# backwards compatibility 
		self.get_a_coef_iterator = self.get_b_coef_iterator = None


	@staticmethod
	def _create_indented_calculation_method(series_functions, a_num_of_coefs, b_num_of_coefs):
		n = var('n')
		l = var('l')
		m = var('m')

		a_coefs_vec = [var('a_{' + str(i) + '}') for i in range(a_num_of_coefs)]
		b_coefs_vec = [var('b_{' + str(i) + '}') for i in range(b_num_of_coefs)]

		a_sym_expr = series_functions[0](n, a_coefs_vec) #lambda x: series_functions[0](x, a_coefs_vec)
		b_sym_expr = series_functions[1](n, b_coefs_vec) #lambda x: series_functions[1](x, b_coefs_vec)

		# We want our expressions to use normal numbers exclusively.
		# to get the rational indents to be normal expressions we'll expand them by a normal factor
		# for any an - bn, one may multiply an and bn as follows:
		# 	an -> an*C
		#	bn -> bn*C*C
		# and the limit L will change to L*C.
		# TODO - add ref to paper
		# We'll calculate a(n + l/m) and b(n + l/m) and find a power of m we can multilfy to an and bn 
		# using this method and get a natural expression.

		indented_an = Poly(a_sym_expr.subs(n, n + l/m), n)
		indented_bn = Poly(b_sym_expr.subs(n, n + l/m), n)

		smallest_m_power_an = 0
		smallest_m_power_bn = 0

		for coeff in indented_an.all_coeffs():
			smallest_m_power_an = min(smallest_m_power_an, coeff.as_powers_dict()[m])
		for coeff in indented_bn.all_coeffs():
			smallest_m_power_bn = min(smallest_m_power_bn, coeff.as_powers_dict()[m])

		inflation_power = abs(min(smallest_m_power_an, smallest_m_power_bn / 2))
		if inflation_power % 1 != 0:
			inflation_power = int(expansion_power) + 1

		inflated_an = indented_an * (m ** inflation_power)
		inflated_bn = indented_bn * (m ** (2 * inflation_power))

		return RationalIndents._get_indented_iterator(inflated_an, a_coefs_vec, l, m), \
			RationalIndents._get_indented_iterator(inflated_bn, b_coefs_vec, l, m)
	
	@staticmethod
	def _get_coeffs(sym_expr, free_vars, coeffs_vec, l, m):
		# assuming sym_expr is poly of n
		coeffs_dict = {}
		for power, coef in sym_expr.as_dict().items():
			# power is a tuple (power,) - unused feature for polys with multiple variables
			power = power[0]

			subs = [(sym, v) for sym, v in zip(coeffs_vec + [l, m], free_vars)]
			coeffs_dict[power] = int(coef.subs(subs))

		return [coeffs_dict[p] if p in coeffs_dict.keys() else 0 for p in range(sym_expr.degree()+1)]

	@staticmethod
	def _get_indented_iterator(sym_expr, coeffs_vec, l, m):
		"""
		To support all of the enumerators, we must supply an iterator function of this form:
			iterator(free_vars, max_runs, start_n=1)
		This function uses the symbolic expression (sympy expression) given, to create such a function. 
		The indentation parameters l, m are appended to the end of free_vars.

		:param sym_expr: A sympy.Poly object for the series
		:param coefs_vec: A list of sympy.var, the free variables used in the expression
		:param l: sympy.var object used for the numerator
		:param m: sympy.var object used for the denominator

		"""

		def indented_iterator(free_vars, max_runs, start_n=1):
			#subs = [(sym, v) for sym, v in zip(coeffs_vec + [l, m], free_vars)]
			coeffs = RationalIndents._get_coeffs(sym_expr, free_vars, coeffs_vec, l, m)
			for i in range(start_n, max_runs):
				val = coeffs[0]
				powers_of_i = 1
				for coef in coeffs[1:]:
					powers_of_i *= i
					val += coef * powers_of_i

				yield val

		return indented_iterator

	def get_num_iterations(self):
		return self.num_iters

	def get_calculation_method(self):
		# TODO - move creation logic to init.
		return RationalIndents._create_indented_calculation_method(
			self.series_functions,
			len(self.seeds[0][0]),
			len(self.seeds[0][1])
			)

	def dump_domain_ranges(self):
		# Not supporting older iterators
		pass

	def get_an_degree(self, coeffs):
		return self.an_deg

	def get_bn_degree(self, coeffs):
		return self.bn_deg

	def get_an_length(self):
		return self.an_length

	def get_bn_length(self):
		return self.bn_length

	def iter_polys(self, primary_looped_domain):
		# primary looped_domain is kind of meaningless here since both domains have
		# the same size
		for seed_an, seed_bn in self.seeds:
			for denom in range(2, self.max_denom+1):
				for numer in range(1, denom):
					if gcd(numer, denom) != 1:
						continue
					# both series are indented with the same number
					yield tuple(seed_an + [numer, denom]), tuple(seed_bn + [numer, denom])

	def split_domains_to_processes(self, number_of_instances):
		# assuming number_of_instances < number of seeds
		seeds_chunks = array_split(self.seeds, number_of_instances)
		sub_domains = []
		for seeds_chunk in seeds_chunks:
			next_instance = deepcopy(self)
			next_instance.seeds = seeds_chunk
			sub_domains.append(next_instance)

		return sub_domains