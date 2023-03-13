import mpmath as mp
from collections import Counter
from functools import reduce
from itertools import chain, combinations_with_replacement
from operator import mul

def get_exponents(degree, order, total_consts):
    if degree == None or order == None:
        raise Exception('degree and order cannot be None')
    return list(c for c in map(Counter, chain.from_iterable(combinations_with_replacement(range(total_consts), i) for i in range(degree + 1)))
                if not any(i for i in c.values() if i > order))

def poly_get(consts, exponents):
    mp.mp.dps = min(c.precision for c in consts)
    values = [mp.mpf(str(c.value)) for c in consts]
    if 1 in values: # solely for backwards-compatibility. We don't need 1 in the DB!
        return None
    return [reduce(mul, (values[i] ** exp[i] for i in range(len(values))), mp.mpf(1)) for exp in exponents]

def poly_eval(poly, coeffs, precisions):
    mp.mp.dps = max(precisions) + 10
    min_prec = min(precisions)
    return min(mp.floor(-mp.log10(abs(mp.fdot(poly, coeffs)))), min_prec)

def poly_verify(consts, degree = None, order = None, relation = None, full_relation = None, exponents = None):
    if full_relation:
        degree = full_relation[0]
        order = full_relation[1]
        relation = full_relation[2:]
    if not exponents:
        exponents = get_exponents(degree, order)
    poly = poly_get(consts, get_exponents(degree, order, len(consts)))
    if not poly:
        return None
    return poly_eval(poly, relation, [c.base.precision for c in consts])

def poly_check(consts, degree = None, order = None, exponents = None):
    if not exponents:
        exponents = get_exponents(degree, order)
    try:
        poly = poly_get(consts, exponents)
        if poly:
            mp.mp.dps = 15 # intentionally low-resolution to quickly try something basic...
            res = mp.pslq(poly)
            if res: # then calculating the substance in what we just found!
                return res, poly_eval(poly, res, [c.base.precision for c in consts]), min_prec
    except ValueError:
        # one of the constants has too small precision, or one constant
        # is small enough that another constant is smaller than its precision.
        # eitherway there's no relation to be found here!
        pass 
    return None, None, None
