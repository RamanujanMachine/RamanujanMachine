from unittest import TestCase
from mobius import GeneralizedContinuedFraction
from mobius import SimpleContinuedFraction
from mobius import find_transform
from mobius import MobiusTransform
from mobius import EfficientGCF
from collections import namedtuple
from series_generators import create_series_from_polynomial, create_series_from_compact_poly
from series_generators import create_series_from_shift_reg
import massey
import mpmath
from sympy import pprint
from sympy import GoldenRatio as phi
from sympy import pi
from sympy import E as e
from sympy import besseli
from sympy import lambdify
from sympy import zeta
import sympy
import data.data
from enumerate_over_gcf import EnumerateOverGCF, LHSHashTable
from enumerate_over_signed_rcf import SignedRcfEnumeration
from convergence_rate import calculate_convergence


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
        compare LHS and RHS.
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
            pprint(sympy.Eq(rhs_sym, rhs.sym_expression(8)))
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
                    rate = calculate_convergence(rhs, lambdify((), t.lhs, 'mpmath')())
                    print("Converged with a rate of {} digits per term".format(mpmath.nstr(rate, 5)))

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
                    rhs = SimpleContinuedFraction.from_irrational_constant(lambdify((), lhs, modules="mpmath"),
                                                                           self.precision // 5)
                    shift_reg = massey.slow_massey(rhs.a_, 199)
                    self.assertLessEqual(len(shift_reg), 20)
                    self.compare(lhs, rhs, self.precision//20)

    def test_alternating_sign_simple_cf(self):
        f_sym = e / (e - 1)
        shift_reg_cmp = [1, 0, -2, 0, 1]
        f_const = lambdify((), f_sym, modules="mpmath")
        with mpmath.workdps(self.precision):
            lhs = f_sym
            rhs = GeneralizedContinuedFraction.from_irrational_constant(f_const, [1, -1]*(self.precision//10))
            self.compare(lhs, rhs, self.precision//20)
            shift_reg = massey.slow_massey(rhs.a_, 199)
            self.assertEqual(len(shift_reg), len(shift_reg_cmp))
            for i in range(len(shift_reg)):
                self.assertEqual(shift_reg[i], shift_reg_cmp[i])

    def test_negative_massey_and_cf(self):
        """
        negative test - there is no polynomial logic behind the pi CF sequence.
        """
        with mpmath.workdps(1000):
            rhs = SimpleContinuedFraction.from_irrational_constant(mpmath.pi, 200)
            shift_reg = massey.slow_massey(rhs.a_, 5657)
            self.assertTrue(len(shift_reg) > 99)

    def known_data_test(self, cf_data):
        """
        test all "known data" that is in data.py
        :param cf_data: a database defined in data.py
        :type cf_data: CFData
        """
        with mpmath.workdps(2000):
            for t in cf_data:
                with self.subTest(test_constant=t):
                    d = cf_data[t]
                    rhs_an = create_series_from_shift_reg(d.rhs_an.shift_reg, d.rhs_an.initials, 400)
                    rhs_bn = create_series_from_shift_reg(d.rhs_bn.shift_reg, d.rhs_bn.initials, 400)
                    rhs = GeneralizedContinuedFraction(rhs_an, rhs_bn)
                    self.compare(d.lhs, rhs, 100)

    def test_known_pi_cf(self):
        """
        test known pi CFs
        """
        self.known_data_test(data.data.pi_cf)

    def test_known_e_cf(self):
        """
        test known e CFs
        """
        self.known_data_test(data.data.e_cf)

    def test_known_zeta_cf(self):
        """
        test known CFs of the zeta function
        """
        self.known_data_test(data.data.zeta_cf)

    def test_weird_stuff(self):
        self.known_data_test(data.data.weird_stuff)

    def test_find_transform(self):
        cases = [(2*pi/3, pi), (1/zeta(3), zeta(3)), ((8*e + 21) / (7 - 5*e), e), (25 + besseli(1, 2), besseli(1, 2))]
        with mpmath.workdps(100):
            for t in cases:
                with self.subTest(find_transform=sympy.pretty(t[0])):
                    y_value = lambdify((), t[0], modules="mpmath")()
                    x_value = lambdify((), t[1], modules="mpmath")()
                    transformation = find_transform(x_value, y_value, 30)
                    sym_transform = transformation.sym_expression(t[1])
                    self.assertEqual(sym_transform, t[0])

    def test_enumeration_over_gcf_hashtable(self):
        hashtable = LHSHashTable(range(3), range(3), mpmath.pi, 1e-7)
        hashtable.save('tmp_test.p')
        hashtable_load = hashtable.load_from('tmp_test.p')
        self.assertEqual(hashtable, hashtable_load)

    def test_enumeration_over_gcf(self):
        enumerator = EnumerateOverGCF(sympy.pi, 4)
        results = enumerator.find_hits([[0, 1, 2]]*2, [[0, 1, 2]]*3, print_results=False)
        self.assertEqual(len(results), 1)
        r = results[0]
        an = create_series_from_compact_poly(r.rhs_an_poly, 1000)
        bn = create_series_from_compact_poly(r.rhs_bn_poly, 1000)
        gcf = GeneralizedContinuedFraction(an, bn)
        t = MobiusTransform(r.lhs_coefs)
        with mpmath.workdps(100):
            lhs_val = mpmath.nstr(gcf.evaluate(), 50)
            rhs_val = mpmath.nstr(t(mpmath.pi), 50)
            self.assertEqual(lhs_val, rhs_val)
            self.assertTrue(t.sym_expression(pi) == 4/pi)

    def test_enumerate_signed_RCF(self):
        enumerator = SignedRcfEnumeration(e, 1, [2, 2], 100, 1)
        with mpmath.workdps(self.precision):
            results = enumerator.find_signed_rcf_conj()
            results, duplicates = enumerator.verify_results(results)
            dups = []
            for key in duplicates.keys():
                for dup in duplicates[key]:
                    dups.append(dup)
            adjusted = [[res[0], res[1], list(res[3])] for res in results+dups]
            self.assertTrue([(e/(e-1)), [1, -1], [1, 0, -2, 0, 1]] in adjusted)

    def test_efficient_gcf(self):
        with mpmath.workdps(100):
            gcf_ref = SimpleContinuedFraction.from_irrational_constant(mpmath.pi, 50)
            val = EfficientGCF(gcf_ref.a_, gcf_ref.b_).evaluate()
            val_ref = gcf_ref.evaluate()
            self.assertEqual(val, val_ref)

    def test_known_new_zeta_values(self):
        """
        test new findings of zeta values in data.py
        """
        with mpmath.workdps(2000):
            for t in data.data.new_zeta_findings:
                with self.subTest(test_constant=t):
                    d = data.data.new_zeta_findings[t]
                    rhs_an = [d.rhs_an(i) for i in range(400)]
                    rhs_bn = [d.rhs_bn(i) for i in range(400)]
                    rhs = GeneralizedContinuedFraction(rhs_an, rhs_bn)
                    self.compare(d.lhs, rhs, 100)
                    rate = calculate_convergence(rhs, lambdify((), d.lhs, 'mpmath')())
                    print("Converged with a rate of {} digits per term".format(mpmath.nstr(rate, 5)))
