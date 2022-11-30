from __future__ import annotations
from collections import namedtuple
from decimal import Decimal, getcontext
from enum import Enum
import mpmath as mp
import numpy as np
from sympy import Poly, Symbol, gcd as sgcd, cancel
from time import time
from typing import List, Tuple
CanonicalForm = Tuple[List[int], List[int]]

n = Symbol('n')


class PCF:
    '''
    a polynomial continued fraction, represented by two Polys a, b:
    a0 + b1 / (a1 + b2 / (a2 + b3 / (...)))
    yes, this is the reverse of wikipedia's convention (i.e. https://en.wikipedia.org/wiki/Generalized_continued_fraction)
    '''
    
    a: Poly
    b: Poly

    def __init__(self: PCF, a_coeffs: List[int], b_coeffs: List[int], auto_deflate: bool = True) -> None:
        '''
        a_coeffs, b_coeffs: lists of integers from the largest power to the smallest power.
        '''
        self.a = Poly(a_coeffs, n)
        self.b = Poly(b_coeffs, n)
        if auto_deflate:
            self.deflate()

    def moving_canonical_form(self: PCF) -> Tuple(Poly, Poly):
        # Should always be real roots. (TODO modify if not!)
        roots = [r for r in self.b.all_roots() if r.is_real]

        # In case b is a constant (has no roots) we still want an to have the canonical roots
        # => (the smallest root to be in (-1,0] )
        # If some of the roots are irrational, it makes the coefficients look ugly, so I decided not to move them.
        # ground_roots is a dict {root:power_of_root} while real_roots is a list of all of the roots including multiplicity
        if len(roots) == 0:
            roots = self.a.real_roots()
            if len(roots) == 0 or len(roots) != sum(self.a.ground_roots().values()):
                return self.a, self.b

        largest_root = max(roots)
        # We want the largest root to be in (-1,0].
        return self.b.compose(Poly(n + largest_root)), self.a.compose(Poly(n + largest_root))

    def inflating_canonical_form(self: PCF) -> Tuple(Poly, Poly):
        top = self.b
        bot = self.a * self.a.compose(Poly(n - 1))
        gcd = sgcd(top, bot)
        return Poly(cancel(top / gcd), n), Poly(cancel(bot / gcd), n)

    def get_canonical_form(self: PCF) -> Tuple(Poly, Poly):
        top, bot = self.inflating_canonical_form()
        return PCF(bot.all_coeffs(), top.all_coeffs()).moving_canonical_form()

    def get_canonical_form_string(self: PCF) -> str:
        a, b = self.get_canonical_form()
        return str(b / a)

    def __str__(self: PCF) -> str:
        return f'a: {self.a.all_coeffs()}\t|\tb: {self.b.all_coeffs()}'

    def is_inflation(self: PCF) -> bool:
        return sgcd(self.b, self.a * self.a.compose(Poly(n - 1))) != 1

    def deflate(self: PCF) -> None:
        deflated: bool = True
        while deflated: # keep going so long as something cancels out
            deflated = False
            a_factors = [factor_tuple[0] for factor_tuple in self.a.factor_list()[1]]
            b_factors = [factor_tuple[0] for factor_tuple in self.b.factor_list()[1]]
            for factor in a_factors:
                if factor in b_factors and factor.compose(Poly(n-1)) in b_factors:
                    self.a = Poly(cancel(self.a / factor))
                    self.b = Poly(cancel(self.b / (factor * factor.compose(Poly(n - 1)))))
                    deflated = True
    
    @staticmethod
    def from_canonical_form(canonical_form: CanonicalForm) -> PCF:
        '''
        Receive the canonical form of a pcf (an := 1 ; bn := bn / (an*a(n+1)))
        and return a pcf of this canonical form.
        Notice there may be many pcfs that fit the same canonical form, this returns just one of them.
        TODO: add link to the doc which explains this
        '''
        n = Symbol('n')
        a = Poly(canonical_form[1], n).compose(Poly(n + 1))
        b = Poly(canonical_form[0], n) * a
        return PCF(a.all_coeffs(), b.all_coeffs())


CALC_JUMP = 256
REDUCE_JUMP = 128
LOG_CALC_JUMP = 7
LOG_REDUCE_JUMP = 6
FR_THRESHOLD = 0.1

class PCFCalc:

    class Convergence(Enum):
        ZERO_DENOM = 0 # now considered an illegal PCF
        NO_FR = 1 # now considered an illegal PCF
        INDETERMINATE_FR = 2
        FR = 3
        RATIONAL = 4
    
    class IllegalPCFException(Exception):
        pass

    class NoFRException(Exception):
        pass
    
    Result = namedtuple('Result', ['value', 'precision', 'last_matrix', 'depth', 'convergence'])

    a: Poly
    b: Poly
    mat: mp.matrix
    depth: int
    
    def __init__(self: PCFCalc, pcf: PCF, prev: List[int] or None = None, depth: int = 0):
        self.a = pcf.a
        self.b = pcf.b
        #self.reduction = 1 # used to exist in the old code, not sure what use we have for this
        self.mat = mp.matrix([prev[0:2], prev[2:4]] if prev else [[self.a(0), 1], [1, 0]])
        self.depth = depth
    
    def reduce(self: PCFCalc):
        gcd = np.gcd.reduce(self.mat)
        #self.reduction *= gcd
        self.mat /= mp.mpf(gcd)
    
    @property
    def value(self: PCFCalc):
        return self.mat[0,0] / self.mat[0,1]
    
    @property
    def precision(self: PCFCalc):
        return mp.floor(-mp.log10(abs(self.value - self.mat[1,0] / self.mat[1,1]))) if all(self.mat[:,1]) else -mp.inf
    
    def check_convergence(self: PCFCalc, fr_list) -> PCFCalc.Convergence:
        p, q = self.mat[0,:]
        val = p / q
        if mp.almosteq(p, 0) or mp.almosteq(val, 0) or mp.pslq([val, 1], tol=mp.power(10, -100)):
            return PCFCalc.Convergence.RATIONAL
        
        if any(abs(fr_list[i + 1] - fr_list[i]) < FR_THRESHOLD for i in range(len(fr_list) - 1)):
            return PCFCalc.Convergence.FR
        
        if any(abs(fr_list[i + 1] - fr_list[i + 2]) > abs(fr_list[i] - fr_list[i + 1]) for i in range(len(fr_list) - 2)):
            return PCFCalc.Convergence.NO_FR
        
        return PCFCalc.Convergence.INDETERMINATE_FR
    
    def run(self: PCFCalc, **kwargs) -> PCFCalc.Result or Exception:
        '''
        Approximate the value of this PCF. Will first calculate to an initial depth,
        and then will procedurally double the depth until the desired precision is obtained.
        
        Accepted kwargs: (others will be ignored)
            'depth': Minimal depth to calculate to, defaults to 10000
            'precision': Minimal precision to obtain, defaults to 50
            'force_fr': Ensure the result has FR if possible (AKA keep going if INDETERMINATE_FR), defaults to True
            'timeout_sec': If nonzero, halt calculation after this many seconds and return whatever you got, no matter what. Defaults to 0
            'timeout_check_freq': Only check for timeout every this many iterations. Defaults to 1000
            'no_exception': Return exceptions (or PCFCalc.Result if possible) instead of raising them (see below). Defaults to False
        
        Exceptions:
            NoFRException: The PCF doesn't converge.
            IllegalPCFException: The PCF has natural roots.
        '''
        # P.S.: The code here is in fact similar to enumerators.FREnumerator.check_for_fr, but
        # the analysis we're doing here is both more delicate (as we allow numerically-indeterminate PCFs),
        # and also less redundant (since we also want the value of the PCF instead of just discarding it for instance)
        def combine(self: PCFCalc, mats, force: bool = False):
            orig = len(mats)
            while len(mats) > 1 and (force or mats[-1][1] >= mats[-2][1]):
                mat1 = mats.pop()
                mat2 = mats.pop()
                mats += [(mat1[0] * mat2[0], mat1[1] + mat2[1])]
            if force or orig - len(mats) > LOG_REDUCE_JUMP:
                gcd = np.gcd.reduce(mats[-1][0])
                mats[-1] = (mats[-1][0] / mp.mpf(gcd), mats[-1][1])
            if force or orig - len(mats) > LOG_CALC_JUMP:
                return mats, mp.log(mp.mpf(np.gcd(*mats[-1][0][0,:]))) / mp.mpf(self.depth) + self.a.degree() * (1 - mp.log(self.depth))
            return mats, # this comma is not a typo! this becomes a 1-tuple
        
        DEFAULTS = {
            'depth': 8192, # at this depth, calculation of one PCF is expected to take about 3 seconds, depending on your machine
            'precision': 50,
            'force_fr': True,
            'timeout_sec': 0,
            'timeout_check_freq': 1000,
            'no_exception': False
        }
        kwargs = {**DEFAULTS, **kwargs}
        getcontext().prec = 2000
        mp.mp.dps = 2000
        fr_list = []
        start = time()
        mats = [(self.mat, self.depth)]
        while self.depth < kwargs['depth']:
            self.depth += 1
            res = combine(self, mats + [(mp.matrix([[self.a(self.depth), self.b(self.depth)], [1, 0]]), 1)])
            if len(res) > 1:
                fr_list += [res[1]]
            mats = res[0]
            if kwargs['timeout_sec'] and self.depth % kwargs['timeout_check_freq'] == 0 and time() - start > kwargs['timeout_sec']:
                break
            if self.depth == kwargs['depth']:
                res = combine(self, mats, True)
                #fr_list += [res[1]]
                self.mat = res[0][0][0]
                
                prec = self.precision # check precision
                if prec == -mp.inf:
                    ex = PCFCalc.IllegalPCFException('continuant denominator zero')
                    if kwargs['no_exception']:
                        return ex
                    raise ex
                if prec < kwargs['precision']:
                    kwargs['depth'] *= 2
                    continue
                
                if kwargs['force_fr']: # check convergence
                    convergence = self.check_convergence(fr_list)
                    if convergence == PCFCalc.Convergence.NO_FR:
                        ex = PCFCalc.NoFRException()
                        if kwargs['no_exception']:
                            return ex
                        raise ex
                    if convergence == PCFCalc.Convergence.INDETERMINATE_FR:
                        kwargs['depth'] *= 2
        
        res = combine(self, mats, True)
        #fr_list += [res[1]]
        self.mat = res[0][0][0]
        value = self.value
        if value and mp.almosteq(0, value):
            value = 0
        
        prec = self.precision
        if prec == -mp.inf:
            ex = PCFCalc.IllegalPCFException('continuant denominator zero')
            if not kwargs['no_exception']:
                raise ex
        
        convergence = self.check_convergence(fr_list)
        if convergence == PCFCalc.Convergence.NO_FR:
            ex = PCFCalc.NoFRException()
            if not kwargs['no_exception']:
                raise ex
        
        return PCFCalc.Result(Decimal(str(value)), 2000 if prec == mp.inf else int(prec), [int(x) for x in self.mat], self.depth, convergence.value)


if __name__ == "__main__":
    b = Poly([-1, 14, -84, 280, -560, 672, -448, 128, 0], n)
    a = Poly([4, -56, 360, -1384, 3476, -5844, 6414, -4170, 1215], n)
    mypcf = PCF(a.all_coeffs(), (b * a.compose(Poly(n + 1))).all_coeffs())
    mypcf.deflate()
    mypcf.moving_canonical_form()
