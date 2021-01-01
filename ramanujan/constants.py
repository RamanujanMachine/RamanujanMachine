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

    mpf_val = mpmath.khinchin   # Hackish trick, used in EnumerationOverGCF.__init__

    def __str__(self):
        return 'K'

    def _latex(self, printer):
        return r"\Kai"


sympy.S.register(Khinchin)
