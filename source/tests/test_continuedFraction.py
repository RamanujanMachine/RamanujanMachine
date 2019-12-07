from unittest import TestCase
from mobius import GeneralizedContinuedFraction
from mobius import SimpleContinuedFraction
from collections import namedtuple
from massey import create_series_from_shift_reg
from massey import create_series_from_polynomial
import massey
import mpmath
from sympy import pprint
# from sympy import GoldenRatio as phi # Commented out due to problems in lambdify.
from sympy import pi
from sympy import E as e
from sympy import besseli
from sympy import lambdify
import sympy
import data.data
phi = (1+sympy.sqrt(5))/2


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

    def compare(self, lhs, rhs, n, print_result=True):
        """

        :param print_result: True will pretty print the comparision
        :param n: number of decimal digits to compare
        :type lhs: sym equation
        :type rhs: GeneralizedContinuedFraction
        """
        rhs_val = mpmath.nstr(rhs.evaluate(), n)
        lhs_f = lambdify((), lhs, modules="mpmath")
        lhs_val = mpmath.nstr(lhs_f(), n)
        if print_result:
            lhs_sym, rhs_sym = sympy.symbols('LHS RHS')
            print("Comparing:")
            pprint(sympy.Eq(lhs_sym, lhs))
            print("With:")
            pprint(sympy.Eq(rhs_sym, rhs.sym_expression(4)))
            if rhs_val == lhs_val:
                print("They Are Equal!\n")
        self.assertEqual(rhs_val, lhs_val)

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
                         lhs=e / (e - 2),
                         rhs_an=create_series_from_polynomial([4, 1], self.depth),
                         rhs_bn=create_series_from_polynomial([-1, -1], self.depth)),
            TestConstant(name='e2',
                         lhs=1 / (e - 2) + 1,
                         rhs_an=create_series_from_polynomial([1], self.depth),
                         rhs_bn=create_series_from_shift_reg([1, 0, -2, 0, 1], [1, -1, 2, -1], self.depth)),
            TestConstant(name='pi1',
                         lhs=4 / (3 * pi - 8),
                         rhs_an=create_series_from_polynomial([3, 3], self.depth),
                         rhs_bn=create_series_from_shift_reg([1, -3, 3, -1], [-1 * 1, -2 * 3, -3 * 5], self.depth)),
            TestConstant(name='phi1',
                         lhs=phi,
                         rhs_an=create_series_from_polynomial([1], self.depth),
                         rhs_bn=create_series_from_polynomial([1], self.depth))
        ]
        with mpmath.workdps(self.precision):
            for t in test_cases:
                with self.subTest(test_constant=t.name):
                    rhs = GeneralizedContinuedFraction(t.rhs_an, t.rhs_bn)
                    self.compare(t.lhs, rhs, self.precision)

    def test_massey_and_creation_of_simple_continued_fractions(self):
        """
        unittests for our regular continued fractions
        """
        rcf_constants = {
            'e': e,
            'bessel_ratio': besseli(1, 2) / besseli(0, 2),
            'phi': phi
        }
        with mpmath.workdps(self.precision):
            for c in rcf_constants:
                with self.subTest(test_constant=c):
                    lhs = rcf_constants[c]
                    rhs = SimpleContinuedFraction(lambdify((), lhs, modules="mpmath"), self.precision // 5)
                    shift_reg = massey.slow_massey(rhs.a_, 199)
                    self.assertLessEqual(len(shift_reg), 20)
                    self.compare(lhs, rhs, self.precision//20)

    def test_negative_massey_and_cf(self):
        """
        negative test - there is no polynomial logic behind the pi CF sequence.
        """
        with mpmath.workdps(1000):
            rhs = SimpleContinuedFraction(mpmath.pi, 200)
            shift_reg = massey.slow_massey(rhs.a_, 5657)
            self.assertTrue(len(shift_reg) > 99)

    def known_data_test(self, cf_data):
        with mpmath.workdps(2000):
            for t in cf_data:
                with self.subTest(test_constant=t):
                    d = cf_data[t]
                    rhs_an = massey.create_series_from_shift_reg(d.rhs_an.shift_reg, d.rhs_an.initials, 400)
                    rhs_bn = massey.create_series_from_shift_reg(d.rhs_bn.shift_reg, d.rhs_bn.initials, 400)
                    rhs = GeneralizedContinuedFraction(rhs_an, rhs_bn)
                    self.compare(d.lhs, rhs, 100)

    def test_known_pi_cf(self):
        self.known_data_test(data.data.pi_cf)

    def test_known_e_cf(self):
        self.known_data_test(data.data.e_cf)

    def test_known_zeta_cf(self):
        self.known_data_test(data.data.zeta_cf)
