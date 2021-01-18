from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficentGCFEnumerator import EfficentGCFEnumerator
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.constants import g_const_dict

# create a LHS table for zeta
saved_hash = 'zeta3.lhs.dept20.db'
lhs_search_limit = 20
lhs = LHSHashTable(
    saved_hash,
    lhs_search_limit,
    [g_const_dict['zeta'](3)])

# define a poly domain to iter through
# Zeta3Domain1 defined an and bn as follows:
# an = (x0 * n + x1) (x2 * n * (n + 1) + x3
# bn = x0 * n^6
# also, the domain has will skip combination we can prove that will not converge.
poly_search_domain = Zeta3Domain1(
    [(2, 2), (1, 1), (-50, 50), (-50, 50)],
    (-10, -1))

# create an enumerator to iter thought the poly domain and compare it to the
# lhs table
enumerator = EfficentGCFEnumerator(
    lhs,
    poly_search_domain,
    [g_const_dict['zeta'](3)]
)

results = enumerator.full_execution()
