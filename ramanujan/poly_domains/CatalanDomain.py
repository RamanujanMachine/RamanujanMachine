from .CartesianProductPolyDomain import CartesianProductPolyDomain
from ..utils.utils import iter_series_items_from_compact_poly


class CatalanDomain(CartesianProductPolyDomain):
	"""
	For catalan conjectures, we're limiting bn to take this form:
		bn = x0 * (2*n+1) ** 4 + b1 * (2*n+1) ** 3

	While an is limited to be a 3rd degree polynomial.
	"""
	def __init__(self, a_coefs_ranges=(0, 0), poly_an_degree=3, b_coefs_ranges=((0, 0), (0, 0)), *args, **kwargs):
		"""
		:param a_coefs_ranges: The range allowed for all an coeffs (the same range for all of them).
		:param b_coefs_ranges: The ranges allowd for bn coeffs. You must specify ranges for both coefs.
		"""
		# b_coef_range are is given blank values. It will be initialized again afterwards
		super().__init__(a_deg=poly_an_degree, b_deg=2, a_coef_range=a_coefs_ranges, b_coef_range=[0, 0], *args, **kwargs)
		self.b_coef_range = b_coefs_ranges

		self._setup_metadata()

	@staticmethod
	def get_calculation_method():
		def bn_iterator(free_vars, max_runs, start_n=0):
			for i in range(start_n, max_runs):
				yield free_vars[0] * (2 * i - 1) ** 4 + free_vars[1] * (2 * i - 1) ** 3

		# an is a standard poly and does not require a special iterator
		return iter_series_items_from_compact_poly, bn_iterator

	def get_bn_degree(self, bn_coefs):
		# bn_coefs is not used since the degree is always 4. Making this calculation immediate and quicker
		return 4
