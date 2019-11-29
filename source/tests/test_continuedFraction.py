from unittest import TestCase
from mobius import GeneralizedContinuedFraction
from mobius import SimpleContinuedFraction
from collections import namedtuple
from massey import create_series_from_shift_reg
from massey import create_series_from_polynomial
import massey
import mpmath


class TestContinuedFracture(TestCase):
    def setUp(self):
        """
            setting up test environment.
            precision is used to for constant generating and comparision
            depth is used for the depth of the continued fraction.
            TODO - is there a way to determine depth <=> precision ratio?
        """
        self.precision = 1000
        self.depth = self.precision * 5  # TODO - find better way to determine depth
        self.e = mpmath.e
        self.pi = mpmath.pi
        self.phi = mpmath.phi

    def test_known_constants(self):
        """
            unittests for our Continued Fraction class.
            all examples were taken from the Ramanujan Machine paper.
            when I didn't know how to describe some series's rules, i used the Massey algorithm to generate it's shift
            register. very useful!
        """
        TestConstant = namedtuple('TestConstant', 'name lhs rhs_an rhs_bn')
        test_cases = [
            TestConstant(name='e1',
                         lhs=lambda: self.e / (self.e - 2),
                         rhs_an=create_series_from_polynomial([4, 1], self.depth),
                         rhs_bn=create_series_from_polynomial([-1, -1], self.depth)),
            TestConstant(name='e2',
                         lhs=lambda: 1 / (self.e - 2) + 1,
                         rhs_an=create_series_from_polynomial([1], self.depth),
                         rhs_bn=create_series_from_shift_reg([1, 0, -2, 0, 1], [1, -1, 2, -1], self.depth)),
            TestConstant(name='pi1',
                         lhs=lambda: 4 / (3 * self.pi - 8),
                         rhs_an=create_series_from_polynomial([3, 3], self.depth),
                         rhs_bn=create_series_from_shift_reg([1, -3, 3, -1], [-1 * 1, -2 * 3, -3 * 5], self.depth)),
            TestConstant(name='phi1',
                         lhs=self.phi,
                         rhs_an=create_series_from_polynomial([1], self.depth),
                         rhs_bn=create_series_from_polynomial([1], self.depth))
        ]
        with mpmath.workdps(self.precision):
            for t in test_cases:
                with self.subTest(test_constant=t.name):
                    rhs = GeneralizedContinuedFraction(t.rhs_an, t.rhs_bn)
                    rhs_val = mpmath.nstr(rhs.evaluate(), n=self.precision)
                    lhs_val = mpmath.nstr(t.lhs(), n=self.precision)
                    # print('comparing:\nlhs = {}\nrhs = {}'.format(lhs_val, rhs_val))      # TODO - optional logging!
                    self.assertEqual(rhs_val, lhs_val)

    def test_massey_and_creation_of_simple_continued_fractions(self):
        """
        unittests for our regular continued fractions
        """
        rcf_constants = {
            'e': self.e,
            'bessel_ratio': lambda: mpmath.besseli(1, 2) / mpmath.besseli(0, 2),
            'phi': self.phi
        }
        with mpmath.workdps(self.precision):
            for c in rcf_constants:
                with self.subTest(test_constant=c):
                    lhs = rcf_constants[c]()
                    rhs = SimpleContinuedFraction(rcf_constants[c], self.precision // 5)
                    shift_reg = massey.slow_massey(rhs.a_, 199)
                    print("\tsimple continued fraction of {}:{}".format(c, rhs))           # TODO - optional logging!
                    print("\tmassey shift register:{}".format(shift_reg, len(shift_reg)))  # TODO - optional logging!
                    self.assertLessEqual(len(shift_reg), 20)
                    lhs_val = mpmath.nstr(lhs, self.precision//20)
                    rhs_val = mpmath.nstr(rhs.evaluate(), self.precision//20)
                    self.assertEqual(lhs_val, rhs_val)

    def test_negative_massey_and_cf(self):
        """
        negative test - there is no polynomial logic behind the pi CF sequence.
        """
        with mpmath.workdps(10000):
            rhs = SimpleContinuedFraction(mpmath.pi, 2000)
            shift_reg = massey.slow_massey(rhs.a_, 5657)
            lhs_val = mpmath.nstr(mpmath.pi(), 1000)
            rhs_val = mpmath.nstr(rhs.evaluate(), 1000)
            # print("\tsimple continued fraction of pi:{}                                   # TODO - optional logging!
            # print("\tmassey shift register len: {}".format(shift_reg, len(shift_reg)))    # TODO - optional logging!
            self.assertEqual(lhs_val, rhs_val)
            self.assertTrue(len(shift_reg) > 999)
