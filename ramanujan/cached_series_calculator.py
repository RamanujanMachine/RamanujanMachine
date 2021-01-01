import collections
from .series_generators import iter_series_items_from_compact_poly

class CachedSeriesCalculator(object):
	"""
	Iterates throught a single series while keeping cache for previus items 
	This is a simpler case for CachedPolySeriesCalculator and is handled in a
	class of its own to reduce redundent actions.
	"""
	def __init__(self, poly_coef):
		self.poly_coef = poly_coef
		self.cached_items = []

	def iter_series_items(self, max_iters=1000):
		"""
		Notes:
		1. if you call this function with start_n and then call it again
		with another start_n, caching will break
		2. start_n
		"""
		yield from self.cached_items

		itered_so_far = len(self.cached_items)
		remaining_iters = max_iters - itered_so_far
		for i in iter_series_items_from_compact_poly(self.poly_coef, max_runs=max_iters, 
			start_n=itered_so_far):
			self.cached_items.append(i)
			yield i

