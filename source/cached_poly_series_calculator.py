import collections
from series_generators import iter_series_items_from_compact_poly

class CachedPolySeriesCalculator(object):
	"""
	Handles iterating through different polynomial series while caching computed items.
	This class meant to address the problem that while calculating RHS combinations 
	every item in the internal nested loop will be computed twice.
	"""
	def __init__(self):
		self.cached_items = {}

	def iter_series_items(self, poly_coef, max_iters=1000):
		if poly_coef in self.cached_items:
			for i in self.cached_items[poly_coef]:
				print('pulled from cache')
				yield i
		else:
			self.cached_items[poly_coef] = []

		remaining_iters = max_iters - len(self.cached_items[poly_coef]) + 1
		for i in iter_series_items_from_compact_poly(poly_coef, max_runs=remaining_iters, 
			start_n=start_n):
			self.cached_items[poly_coef].append(i)
			yield i