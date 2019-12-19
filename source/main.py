import mobius
import massey
import mpmath
import time
import sympy
import itertools

def find_transform_slow(x, y, limit):
    def sort_l(l):
        return sum([abs(s) for s in l])

    vals = [i for i in range(-limit, limit)]
    options = [vals, vals, vals, vals]
    values = [e for e in itertools.product(*options)]
    values_sorted = sorted(values, key=sort_l)
    iters = 0
    for v in values_sorted:
        iters += 1
        try:
            res = (v[0] * x + v[1]) / (v[2] * x + v[3])
        except ZeroDivisionError:
            continue
        if abs(res - y) < 1e-10:
            return v, iters
    return [], iters

def f1():
    with mpmath.workdps(100):
        x = mpmath.pi()
        an = massey.create_series_from_shift_reg([1, -2, 1], [-1, -4], 200)
        bn = massey.create_series_from_shift_reg([1, -3, 3, -1], [-2, -7, -9], 200)
        #an = massey.create_series_from_shift_reg([1, -2, 1], [-4, -7], 200)
        #bn = massey.create_series_from_shift_reg([1, -3, 3, -1], [-9, -20, -35], 200)
        gcf = mobius.GeneralizedContinuedFraction(an, bn)
        # gcf.print(5)
        y = gcf.evaluate()
        t = mobius.find_transform(x, y, 20)
        if t is not None:
            t.pprint()
        else:
            print("no solution found")


def check_const(f):
    with mpmath.workdps(20000):
        cf = mobius.SimpleContinuedFraction.from_irrational_constant(f, 200)
        cf.print(8)
        massey.massey_check(cf.a_)

f1()
#res = mobius.find_transform(mpmath.pi, (mpmath.pi*2)/3, 20)
#if res is not None:
#    print('Found a Solution:')
#    res.pprint(sympy.pi)
#pi = sympy.pi
#e = sympy.E
#f_sym = ((pi + 1) / 2) + (pi / (e**(2*pi)-1))
#sympy.pprint(f_sym)
#check_const(sympy.lambdify((), f_sym, modules="mpmath"))