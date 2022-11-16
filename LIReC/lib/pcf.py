from __future__ import annotations
from sympy import Poly, Symbol, gcd as sgcd, cancel
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
        """
        a_coeffs, b_coeffs: lists of integers from the largest power to the smallest power.
        """
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
        """
        Receive the canonical form of a pcf (an := 1 ; bn := bn / (an*a(n+1)))
        and return a pcf of this canonical form.
        Notice there may be many pcfs that fit the same canonical form, this returns just one of them.
        TODO: add link to the doc which explains this
        """
        n = Symbol('n')
        a = Poly(canonical_form[1], n).compose(Poly(n + 1))
        b = Poly(canonical_form[0], n) * a
        return PCF(a.all_coeffs(), b.all_coeffs())


if __name__ == "__main__":
    b = Poly([-1, 14, -84, 280, -560, 672, -448, 128, 0], n)
    a = Poly([4, -56, 360, -1384, 3476, -5844, 6414, -4170, 1215], n)
    mypcf = PCF(a.all_coeffs(), (b * a.compose(Poly(n + 1))).all_coeffs())
    mypcf.deflate()
    mypcf.moving_canonical_form()
