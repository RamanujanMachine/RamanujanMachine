from functools import reduce
from sympy.polys.polytools import poly_from_expr
import os
from db.config import configuration
from db.jobs.job_poly_pslq import execute_job, transpose
from db.lib.ramanujan_db import RamanujanDB
from db.lib.models import PcfCanonicalConstant
from db.lib.pcf import *

PRECISION_LIMIT = 20

def main():
    os.makedirs(os.path.join(os.getcwd(), 'logs'), exist_ok=True)
    config = configuration['analyze_pcfs']
    db_handle = RamanujanDB()
    pcfs = (PCF(*(poly_from_expr(poly)[0].all_coeffs() for poly in pcf)) for pcf in config['pcfs'])
    successful, unsuccessful = db_handle.add_pcfs(pcfs)#[],{'Already exist':pcfs}#
    unsuccessful['Too imprecise'] = []
    for key in unsuccessful:
        if key == 'Already exist': # still want to test these...
            canonical_forms = [[poly.all_coeffs() for poly in pcf.get_canonical_form()] for pcf in unsuccessful[key]]
            for pcf in db_handle.session.query(PcfCanonicalConstant):
                if [pcf.P, pcf.Q] in canonical_forms: # TODO theoretically this can be implemented as a filter but i can't figure it out right now...
                    if pcf.base.precision > PRECISION_LIMIT:
                        successful += [pcf]
                    else:
                        unsuccessful['Too imprecise'] += [pcf]
        elif unsuccessful[key]:
            polys = [(pcf.P, pcf.Q) for pcf in unsuccessful[key]]
            pcfs_string = reduce(lambda a,b: a+'\n\t'+str(b), polys[1:], str(polys[0]))
            print(f'Could not add {len(unsuccessful[key])} pcfs with reason "{key}":\n\t{pcfs_string}\n')
    db_handle.session.close()
    if 'PcfCanonical' in config['subdivide'] and not successful:
        print('None of the configued PCFs can be tested! Aborting...')
    else: # bulk doesn't matter, it just exists due to the limitation of how EXECUTE_NEEDS_ARGS works
        print(f'Testing relations with {len(successful)} PCFs')
        execute_job(transpose([successful]), config.get('subdivide', None), config.get('degree', None), debug_pslq = config.get('debug_pslq', False), manual = True)

if __name__ == '__main__':
    '''
    Example config:
    configuration = {
        ...
        'analyze_pcfs': {
            'degree': (2, 1), 'debug_pslq': False, 'subdivide': {
                'PcfCanonical': { 'count': 1, 'balanced_only': True },
                'Named': { 'count': 1 }
            },
            'pcfs': [
                ('2*n**5 + 7*n**4 + 14*n**3 + 16*n**2 + 9*n + 2', '-n**10 - 2*n**9 - n**8'),
                ('2*n**5 + 7*n**4 + 19*n**3 + 25*n**2 + 16*n + 4', '-n**10 - 2*n**9 - n**8'),
                ('2*n**5 + 7*n**4 + 19*n**3 + 30*n**2 + 21*n + 5', '-n**10 - 2*n**9 - n**8'),
                ('2*n**5 + 8*n**4 + 12*n**3 + 10*n**2 + 5*n + 1', '-n**10 - 3*n**9 - 2*n**8'),
                ('2*n**5 + 8*n**4 + 14*n**3 + 16*n**2 + 9*n + 2', '-n**10 - 3*n**9 - 2*n**8'),
                ('2*n**5 + 8*n**4 + 18*n**3 + 22*n**2 + 13*n + 3', '-n**10 - 3*n**9 - 2*n**8'),
                ('2*n**5 + 8*n**4 + 18*n**3 + 28*n**2 + 19*n + 6', '-n**10 - 3*n**9 - 2*n**8'),
                ('2*n**5 + 8*n**4 + 24*n**3 + 34*n**2 + 23*n + 6', '-n**10 - 3*n**9 - 2*n**8'),
                ('2*n**5 + 10*n**4 + 14*n**3 + 16*n**2 + 9*n + 2', '-n**10 - 5*n**9 - 4*n**8'),
                ('2*n**5 + 10*n**4 + 26*n**3 + 34*n**2 + 21*n + 5', '-n**10 - 5*n**9 - 4*n**8'),
                ('2*n**5 + 11*n**4 + 34*n**3 + 46*n**2 + 29*n + 7', '-n**10 - 6*n**9'),
                ('2*n**5 + 13*n**4 + 10*n**3 + 10*n**2 + 5*n + 1', '-n**10 - 8*n**9'),
                ('2*n**3 + 3*n**2 + 3*n + 1', '-n**6'),
                ('2*n**3 + 3*n**2 + 11*n + 5', '-n**6'),
                ('2*n**3 + 3*n**2 + 27*n + 13', '-n**6'),
                ('2*n**3 + 3*n**2 + 51*n + 25', '-n**6'),
                ('2*n**3 + 3*n**2 + 53*n + 26', '-n**6'),
                ('2*n**3 + 3*n**2 + 81*n + 40', '-n**6'),
                ('2*n**3 + 3*n**2 + 83*n + 41', '-n**6'),
                ('2*n**3 + 3*n**2 + 123*n + 61', '-n**6'),
                ('2*n**3 + 3*n**2 + 171*n + 85', '-n**6'),
                ('2*n**3 + 3*n**2 + 227*n + 113', '-n**6'),
                ('3*n**3 + 3*n**2 + 3*n + 1', '-2*n**6'),
                ('3*n**3 + 6*n**2 + 6*n + 2', '-2*n**6'),
                ('4*n**3 + 3*n**2 + 3*n + 1', '-3*n**6'),
                ('4*n**3 + 9*n**2 + 9*n + 3', '-3*n**6'),
                ('4*n**3 + 6*n**2 + 4*n + 1', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 6*n + 2', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 12*n + 5', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 22*n + 10', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 32*n + 15', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 36*n + 17', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 54*n + 26', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 76*n + 37', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 102*n + 50', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 106*n + 52', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 132*n + 65', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 162*n + 80', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 166*n + 82', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 204*n + 101', '-4*n**6'),
                ('4*n**3 + 6*n**2 + 246*n + 122', '-4*n**6'),
                ('5*n**3 + 3*n**2 + 3*n + 1', '-4*n**6'),
                ('5*n**3 + 12*n**2 + 12*n + 4', '-4*n**6'),
                ('5*n**3 + 9*n**2 + 9*n + 3', '-6*n**6'),
                ('6*n**3 + 9*n**2 + 5*n + 1', '-n**6'),
                ('6*n**3 + 3*n**2 + 3*n + 1', '-5*n**6'),
                ('6*n**3 + 15*n**2 + 15*n + 5', '-5*n**6'),
                ('6*n**3 + 6*n**2 + 6*n + 2', '-8*n**6'),
                ('6*n**3 + 12*n**2 + 12*n + 4', '-8*n**6'),
                ('6*n**3 + 9*n**2 + 9*n + 3', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 23*n + 10', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 33*n + 15', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 81*n + 39', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 115*n + 56', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 125*n + 61', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 153*n + 75', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 159*n + 78', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 233*n + 115', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 243*n + 120', '-9*n**6'),
                ('6*n**3 + 9*n**2 + 249*n + 123', '-9*n**6'),
                ('7*n**3 + 3*n**2 + 3*n + 1', '-6*n**6'),
                ('7*n**3 + 18*n**2 + 18*n + 6', '-6*n**6'),
                ('8*n**3 + 3*n**2 + 3*n + 1', '-7*n**6'),
                ('8*n**3 + 21*n**2 + 21*n + 7', '-7*n**6'),
                ('9*n**3 + 3*n**2 + 3*n + 1', '-8*n**6'),
                ('9*n**3 + 24*n**2 + 24*n + 8', '-8*n**6'),
                ('9*n**3 + 31*n**2 + 34*n + 12', '-8*n**6'),
                ('10*n**3 + 3*n**2 + 3*n + 1', '-9*n**6'),
                ('10*n**3 + 27*n**2 + 27*n + 9', '-9*n**6'),
                ('34*n**3 + 51*n**2 + 27*n + 5', '-n**6'),
                ('3*n**4 + 4*n**3 + 6*n**2 + 4*n + 1', '-2*n**8'),
                ('3*n**4 + 8*n**3 + 12*n**2 + 8*n + 2', '-2*n**8'),
                ('4*n**4 + 4*n**3 + 6*n**2 + 4*n + 1', '-3*n**8'),
                ('4*n**4 + 12*n**3 + 18*n**2 + 12*n + 3', '-3*n**8'),
                ('4*n**4 + 8*n**3 + 12*n**2 + 8*n + 2', '-4*n**8'),
                ('4*n**4 + 8*n**3 + 20*n**2 + 16*n + 6', '-4*n**8'),
                ('5*n**4 + 4*n**3 + 6*n**2 + 4*n + 1', '-4*n**8'),
                ('5*n**4 + 16*n**3 + 24*n**2 + 16*n + 4', '-4*n**8'),
                ('10*n**4 + 4*n**3 + 6*n**2 + 4*n + 1', '-9*n**8'),
                ('10*n**4 + 36*n**3 + 54*n**2 + 36*n + 9', '-9*n**8'),
                ('2*n**4 + 4*n**3 + 6*n**2 + 4*n + 1', '-n**8'),
                ('2*n**4 + 4*n**3 + 10*n**2 + 8*n + 3', '-n**8'),
                ('4*n**4 + 8*n**3 + 12*n**2 + 8*n + 2', '-4*n**8'),
                ('4*n**4 + 8*n**3 + 20*n**2 + 16*n + 6', '-4*n**8'),
                ('6*n**4 + 12*n**3 + 18*n**2 + 12*n + 3', '-9*n**8'),
                ('6*n**4 + 12*n**3 + 30*n**2 + 24*n + 9', '-9*n**8'),
                ('8*n**4 + 16*n**3 + 24*n**2 + 16*n + 4', '-16*n**8'),
                ('8*n**4 + 16*n**3 + 40*n**2 + 32*n + 12', '-16*n**8'),
                ('16*n**5 + 40*n**4 + 50*n**3 + 35*n**2 + 13*n + 2', '-64*n**10'),
                ('16*n**5 + 40*n**4 + 56*n**3 + 44*n**2 + 18*n + 3', '-64*n**10'),
                ('2*n**5 + 5*n**4 + 10*n**3 + 10*n**2 + 5*n + 1', '-n**10'),
                ('2*n**5 + 5*n**4 + 22*n**3 + 28*n**2 + 15*n + 3', '-n**10'),
                ('2*n**5 + 5*n**4 + 42*n**3 + 58*n**2 + 45*n + 13', '-n**10'),
                ('16*n**5 + 40*n**4 + 120*n**3 + 140*n**2 + 90*n + 23', '-64*n**10'),
                ('-16*n**5 - 40*n**4 - 120*n**3 - 140*n**2 - 90*n - 23', '-64*n**10'),
                ('-2*n**5 - 5*n**4 - 42*n**3 - 58*n**2 - 45*n - 13', '-n**10'),
                ('16*n**5 + 40*n**4 + 50*n**3 + 35*n**2 + 13*n + 2', '-64*n**10'),
                ('16*n**5 + 40*n**4 + 56*n**3 + 44*n**2 + 18*n + 3', '-64*n**10'),
                ('-2*n**5 - 5*n**4 - 22*n**3 - 28*n**2 - 23*n - 7', '-n**10'),
                ('-2*n**5 - 5*n**4 - 10*n**3 - 10*n**2 - 5*n - 1', '-n**10'),
                ('2*n**5 + 5*n**4 + 22*n**3 + 28*n**2 + 15*n + 3', '-n**10'),
                ('-16*n**5 - 40*n**4 - 56*n**3 - 44*n**2 - 18*n - 3', '-64*n**10'),
                ('-16*n**5 - 40*n**4 - 50*n**3 - 35*n**2 - 13*n - 2', '-64*n**10'),
                ('2*n**5 + 5*n**4 + 42*n**3 + 58*n**2 + 45*n + 13', '-n**10'),
                ('16*n**5 + 40*n**4 + 120*n**3 + 140*n**2 + 90*n + 23', '-64*n**10'),
                ('2*n**5 + 5*n**4 + 10*n**3 + 10*n**2 + 5*n + 1', '-n**10'),
                ('2*n**5 + 5*n**4 + 22*n**3 + 28*n**2 + 15*n + 3', '-n**10'),
                ('2*n**5 + 5*n**4 + 22*n**3 + 28*n**2 + 23*n + 7', '-n**10'),
                ('2*n**5 + 5*n**4 + 42*n**3 + 58*n**2 + 45*n + 13', '-n**10')
            ]
        }
        ...
    }
    '''
    main()
