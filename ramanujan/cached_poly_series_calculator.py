import collections
from .series_generators import iter_series_items_from_compact_poly

class CachedPolySeriesCalculator(object):
	"""
	Handles iterating through different polynomial series while caching computed items.
	This class meant to address the problem that while calculating RHS combinations 
	every item in the internal nested loop will be computed twice.
	"""
	def __init__(self):
		self.cached_items = {}

	def iter_series_items(self, poly_coef, max_iters=1000, series_iterator=iter_series_items_from_compact_poly):
		"""
		Notes:
		1. if you call this function with start_n and then call it again
		with another start_n, caching will break
		2. start_n
		"""
		cache = []

		if poly_coef in self.cached_items:
			cache = self.cached_items[poly_coef]
			yield from cache[:max_iters]
		else:
			self.cached_items[poly_coef] = []

		itered_so_far = len(cache)
		remaining_iters = max_iters - itered_so_far
		
		for i in series_iterator(poly_coef, max_runs=max_iters, 
			start_n=itered_so_far):
			self.cached_items[poly_coef].append(i)
			yield i

	def iter_family(self, coef_iter, max_iters=1000, series_iterator=iter_series_items_from_compact_poly):
#		if self.cached_items == {}:
		# nothing in cache so far
		for coef in coef_iter:
			yield coef, self.iter_series_items(coef, max_iters, series_iterator=series_iterator)
		# else:
		# 	for coef in self.cached_items:
		# 		yield coef, self.iter_series_items(coef, max_iters, series_iterator=iter_series_items_from_compact_poly)