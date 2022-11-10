from logging import getLogger
from logging.config import fileConfig
import numpy as np
from functools import reduce
from operator import mul
from sympy import Poly, Symbol
from sympy.core.numbers import Integer, Rational
from LIReC.lib.pcf import PCF
from LIReC.lib.db_access import LIReC_DB
from LIReC.config import configuration

fileConfig('logging.config', defaults={'log_filename': 'generate_random'})

LOGGER_NAME = 'job_logger'
POLY_MAX_DEGREE = 3
POLY_MAX_COEFF = 50
POLY_MAX_ABS_ROOT = 5
BULK_SIZE = 5000

def get_degrees(max_deg, num_denom_factor):
    if num_denom_factor == None:
        return np.random.randint(max_deg + 1), np.random.randint(max_deg + 1)
    
    factor, strict = num_denom_factor
    low_deg = np.random.randint(max_deg + 1)
    if strict:
        high_deg = low_deg * abs(factor)
    else:
        high_deg = np.random.randint(abs(factor) * low_deg, abs(factor) * max_deg + 1)

    return (low_deg, high_deg) if factor > 0 else (high_deg, low_deg)

def has_natural_roots(p: Poly):
    return any(r for r in p.real_roots() if isinstance(r, Integer) and r > 0)

def generate_pcf(max_deg, max_coeff, max_abs_root, num_denom_factor):
    coeffs = range(-max_coeff, max_coeff)
    coeffs_nonzero = [x for x in coeffs if x != 0]
    a_deg, b_deg = get_degrees(max_deg, num_denom_factor)
    a, b = None, None
    n = Symbol('n')
    natural = True # no natural roots! if a has natural roots the PCF diverges! (i think)
    while natural: # reminder that np.concatenate requires the arrays to be inside a tuple
        a = np.concatenate((np.random.choice(coeffs_nonzero, 1), np.random.choice(coeffs, a_deg))).tolist()
        natural = has_natural_roots(Poly(a, n))
    
    coeffs = range(-max_abs_root, max_abs_root)
    coeffs_nonzero = [x for x in coeffs if x != 0]
    natural = True # no natural roots! if b has natural roots then the PCF is trivially rational, thus not interesting
    while natural: # b must have only rational roots. idk why we choose to limit it like this but i don't question it
        b_1 = np.random.choice(coeffs_nonzero, b_deg + 1)
        b_2 = np.random.choice(coeffs, b_deg)
        b = b_1[b_deg] * reduce(mul, [Poly([b_1[i], b_2[i]], n) for i in range(b_deg)], Poly(1, n))
        natural = has_natural_roots(b)
    
    b = b.all_coeffs()
    print((a,b))
    return PCF(a, b)

def execute_job(bulk=0, max_deg=-1, max_coeff=-1, max_abs_root=-1, num_denom_factor=None):
    # yes this was originally a job, but since this requires superuser permissions it got... "promoted"?
    bulk = bulk if bulk else BULK_SIZE
    max_deg = max_deg if max_deg > 0 else POLY_MAX_DEGREE
    max_coeff = max_coeff if max_coeff > 0 else POLY_MAX_COEFF
    max_abs_root = max_abs_root if max_abs_root > 0 else POLY_MAX_ABS_ROOT

    getLogger(LOGGER_NAME).info(f'starting to generate cfs randomly: {bulk}, {max_deg}, {max_coeff}, {num_denom_factor}')
    db = LIReC_DB()
    _, unsuccessful = db.add_pcfs(generate_pcf(max_deg, max_coeff, max_abs_root, num_denom_factor) for i in range(bulk))
    print({k : len(unsuccessful[k]) for k in unsuccessful}) # should never get illegal pcf's here
    getLogger(LOGGER_NAME).info('finished generating')
    db.session.close()

def main():
    '''
    Example config:
    configuration = {
        ...
        'add_random_pcfs': {
            'bulk': 100, 'max_deg': 3, 'max_coeff': 50, 'max_abs_root': 3, 'num_denom_factor': (2, True)
        },
        ...
    }
    '''
    if 'add_random_pcfs' in configuration:
        execute_job(**configuration['add_random_pcfs'])
    else:
        execute_job() # with defaults...

if __name__ == '__main__':
    main()
