from CartesianProductPolyDomain import * 
from itertools import product 

class FixedBNPolyDomain(CartesianProductPolyDomain):
	def __init__(self, a_deg, a_coef_range, b_coefs, force_lead_a_pos=True,
		*args, **kwargs):
		self.a_deg = a_deg
		self.a_coef_range = a_coef_range
		self.b_coefs = b_coefs
		self.b_deg = len(b_coefs) - 1
		self.force_lead_a_pos = force_lead_a_pos

		self.an_length = self.get_an_length()
		self.bn_length = 1 # Duh
		self.num_iterations = self.an_length * self.bn_length

		self.an_domain_range, self.bn_domain_range = self.dump_domain_ranges()

	def expand_var_ranges_to_domain(self, var_ranges):
		return [ \
			[i for i in range(coef[0], coef[1] + 1)] \
				for classoef in var_ranges
		]

	def dump_domain_ranges(self):
		# we might not want the first coef to be 0
		an_domain = [[i for i in range(self.a_coef_range[0], self.a_coef_range[1] + 1)]]
		if self.force_lead_a_pos and 0 in an_domain[0]:
			an_domain[0] = an_domain[0][an_domain[0].index(0) + 1:]


		for i in range(self.a_deg):
			an_domain.append(
				[i for i in range(self.a_coef_range[0], self.a_coef_range[1]+1)])

		return an_domain, [[i] for i in self.b_coefs]