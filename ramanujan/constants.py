import sympy
from sympy.core.compatibility import with_metaclass
from sympy.core.singleton import Singleton
from sympy.core import NumberSymbol
import mpmath


class Khinchin(with_metaclass(Singleton, NumberSymbol)):
    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = None
    is_number = True

    mpf_val = mpmath.khinchin  # Hackish trick, used in EnumerationOverGCF.__init__

    def __str__(self):
        return "K"

    def _latex(self, printer):
        return r"\Kai"


sympy.S.register(Khinchin)

g_N_verify_terms = (
    1000  # number of CF terms to calculate in __refine_results. (verify hits)
)
g_N_verify_compare_length = (
    100  # number of digits to compare in __refine_results. (verify hits)
)
g_N_verify_dps = 2000  # working decimal precision in __refine_results. (verify hits)
g_N_initial_search_terms = (
    32  # number of CF terms to calculate in __first_enumeration (initial search)
)
g_N_initial_key_length = (
    10  # number of digits to compare in __first_enumeration (initial search)
)
g_N_initial_search_dps = (
    50  # working decimal precision in __refine_results. (verify hits)
)

# math constants:
g_const_dict = {
    "zeta": sympy.zeta,
    "e": sympy.E,
    "pi": sympy.pi,
    "pi_sqared": sympy.pi**2,
    "catalan": sympy.Catalan,
    "golden_ratio": sympy.GoldenRatio,
    "khinchin": sympy.S.Khinchin,
    "euler-mascheroni": sympy.EulerGamma,
    "pi-acosh_2": sympy.pi * sympy.acosh(2),
    "polygamma": sympy.polygamma,
}
