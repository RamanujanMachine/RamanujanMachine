class CachedPolySeriesCalculator(object):
	"""
	Handles iterating through different polynomial series while caching computed items.
	This class meant to address the problem that while calculating RHS combinations 
	every item in the internal nested loop will be computed twice.
	"""
	pass