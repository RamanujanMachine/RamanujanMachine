from abc import ABCMeta

class AbstractPolyDomains(metaclass=ABCMeta):
	"""
	This is an abstract class, meant to represent the polynomial domains
	iterated through a run.

	A poly domain is defined here as two families of polynomials. Each family
	is defined by a function with variables, referred to as coefs. Each
	coef has it allowed range size.
	This class will generate all coefs under those restrictions. Each decanted 
	class must also create a calculation method to use those coefs

	Examples:
		A(n) = a*n^2 + b*n + c
		B(n) = d*n^2 + e*n + f
		is classical a poly domain, where both a,b has 6 deg of freedom 

		A(n) = a*n^2 + b*n + c
		B(n) = A(n) * A(n-1)
		This combination only has 3 degrees of freedom. Also, given a,b,c
		the calculation for A and B are different. The way to calculate each
		series by given constants, is the calculation method
	"""
	def __init__(self):
		pass

	def iter_polys(self, primary_looped_domain):
		'''
			On other locations, one might want to choose the nested
			oreder between both poylnomials (for caching the smaller one)
			primary_looped_domain is ment to allow a user to choose the nesting
			order

			primary_looped_domain accepts 'an' or 'bn'
		'''
		pass

	def get_num_iterations(self):
		pass

	def get_calculation_method(self):
		pass

	def dump_domain_ranges(self):
		"""
			Backwards compatibility - some enumerators except this format
		"""
		pass