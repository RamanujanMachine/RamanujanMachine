from unittest import TestCase
from mobius import ContinuedFraction
from collections import namedtuple
import gen_consts
from massey import create_series_from_shift_reg
from massey import create_series_from_polynomial
import mpmath

TestConstant = namedtuple('TestConstant', 'name lhs rhs_an rhs_bn')


class TestContinuedFracture(TestCase):
    def setUp(self):
        """
            setting up test environment.
            precision is used to for constant generating and comparision
            depth is used for the depth of the continued fraction.
            TODO - is there a way to determine depth <=> precision ratio?
        """
        self.precision = 30
        self.depth = 100
        self.cases = {}

        self.e = gen_consts.gen_e_const(self.precision)
        self.pi = gen_consts.gen_pi_const(self.precision)
        self.phi = gen_consts.gen_phi_const(self.precision)

    def test_known_constants(self):
        """
            unittests for our Continued Fraction class.
            all examples were taken from the Ramanujan Machine paper.
            when I didn't know how to describe some series's rules, i used the Massey algorithm to generate it's shift
            register. very useful!
        """
        test_cases = [
            TestConstant(name='e1',
                         lhs=self.e / (self.e - 2),
                         rhs_an=create_series_from_polynomial([4, 1], self.depth),
                         rhs_bn=create_series_from_polynomial([-1, -1], self.depth)),
            TestConstant(name='e2',
                         lhs=1 / (self.e - 2) + 1,
                         rhs_an=create_series_from_polynomial([1], self.depth),
                         rhs_bn=create_series_from_shift_reg([1, 0, -2, 0, 1], [1, -1, 2, -1], self.depth)),
            TestConstant(name='pi1',
                         lhs=4 / (3 * self.pi - 8),
                         rhs_an=create_series_from_polynomial([3, 3], self.depth),
                         rhs_bn=create_series_from_shift_reg([1, -3, 3, -1], [-1 * 1, -2 * 3, -3 * 5], self.depth)),
            TestConstant(name='phi1',
                         lhs=self.phi,
                         rhs_an=create_series_from_polynomial([1], self.depth),
                         rhs_bn=create_series_from_polynomial([1], self.depth))
        ]
        for t in test_cases:
            with self.subTest(test_constant=t.name):
                rhs = ContinuedFraction(t.rhs_an, t.rhs_bn)
                rhs_val = mpmath.nstr(rhs.evaluate(), n=self.precision)
                lhs_val = mpmath.nstr(t.lhs, n=self.precision)
                # print('comparing:\nlhs = {}\nrhs = {}'.format(lhs_val, rhs_val))
                self.assertEqual(rhs_val, lhs_val)
