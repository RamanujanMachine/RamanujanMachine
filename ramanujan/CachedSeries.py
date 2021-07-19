from ramanujan.utils.utils import iter_series_items_from_compact_poly


class CachedSeries(object):
	"""
	Handles iterating through polynomial series while caching computed items.
	"""
	def __init__(self, poly_coefs, series_iterator=iter_series_items_from_compact_poly):
		self.poly_coefs = poly_coefs
		self.series_iterator = series_iterator
		self.cache = []

	def iter_series_items(self, max_iters=1000):
		yield from self.cache[:max_iters]
		
		# Calculate new items only if needed
		itered_so_far = len(self.cache)
		for i in self.series_iterator(self.poly_coefs, max_runs=max_iters, start_n=itered_so_far):
			self.cache.append(i)
			yield i
