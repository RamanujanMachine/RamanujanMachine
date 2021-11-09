from abc import ABCMeta


class AbstractPolyDomains(metaclass=ABCMeta):
	"""
	This is an abstract class, meant to represent the polynomial domains iterated through a run.

	A poly domain is defined here as two families of polynomials. Each family is defined by a function with variables,
	referred to as coefficients Each coefficient has it allowed range size.
	This class will generate all coefficients under those restrictions.

	Each decanted class must create a calculation method to use those coefficients

	Examples:
		A(n) = a*n^2 + b*n + c
		B(n) = d*n^2 + e*n + f
		is classical a poly domain, where both a,b has 3 degrees of freedom, giving us 6 degrees of freedom in total.

		A(n) = a*n^2 + b*n + c
		B(n) = A(n) * A(n-1)
		This combination only has 3 degrees of freedom. Also, given a,b,c
		the calculation for A and B are different. The way to calculate each
		series by given constants, is the calculation method
	"""
	def __init__(self):
		pass

	def iter_polys(self, primary_looped_domain):
		"""
		On other locations, one might want to choose the nested order between both polynomials
		(for caching the smaller one)
		primary_looped_domain is meant to allow a user to choose the nesting order

		primary_looped_domain accepts 'an' or 'bn'
		"""
		pass

	def get_num_iterations(self):
		pass

	@staticmethod
	def get_calculation_method():
		pass

	def dump_domain_ranges(self):
		"""
		Backwards compatibility - some enumerators except this format.
		returns a list that contains every possible value for each coefficient.
		so [[1,2,3],[1,2,3],[1,2,3]]
		means that each coefficient may accept the values 1, 2 or 3
		"""
		pass
