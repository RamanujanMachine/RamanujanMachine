import mobius
import massey
import mpmath
import time


def find_cn(b_):
    c_ = [1]
    for i in range(1, len(b_)):
        c_.append(1/(b_[i] * c_[i-1]))
    return c_


def create_simple_from_generalized(a_, b_):
    c_ = find_cn(b_)
    a_new_ = [a_[0]]
    for i in range(1, len(a_)):
        a_new_.append(a_[i] * c_[i])
    return a_new_


with mpmath.workdps(10000):
    f = lambda: mpmath.e / (mpmath.e - 2)
    CF = mobius.GeneralizedContinuedFraction(massey.create_series_from_polynomial([4, 1], 200),
                                             massey.create_series_from_polynomial([-1, -1], 200))
    RCF = mobius.SimpleContinuedFraction(f, 2000)
    #print(mpmath.nstr(CF.evaluate(), 50))
    #print(mpmath.nstr(RCF.evaluate(20000), 50))
    print(mpmath.nstr(RCF.evaluate(0), 50))
    #print(CF.mobius)
    #print(RCF.mobius)
    print(RCF.a_)
    #print(CF.a_)
    #print(CF.b_)
    #print(create_simple_from_generalized(CF.a_, CF.b_))
    assert mpmath.nstr(CF.evaluate(), 100) == mpmath.nstr(f(), 100)
    massey.massey_check(RCF.a_)
    print(CF.a_)
    massey.massey_check(CF.a_)
    print(CF.b_)
    massey.massey_check(CF.b_)

    a_ = []

