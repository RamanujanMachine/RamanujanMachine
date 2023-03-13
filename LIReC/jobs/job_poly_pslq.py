'''
Finds polynomial relations between constants in LIReC, using PSLQ.

Configured as such:
'degree' + 'order':
    Two integers. All relations are structured like multivariate polynomials over the constants and CFs,
    of degree 'degree' with a maximum exponent of 'order'. For example, a 2-variable polynomial of
    degree 2 and order 1 will be of the form a+bx+cy+dxy (note the lack of x^2 and y^2), and
    a 4-variable polynomial of degree (3,1) will be of the form:
        a + bx+cy+dz+ew + fxy+gxz+hxw+iyz+jyw+kzw + lxyz+mxyw+nxzw+oyzw
    Note here the lack of any single variable with an exponent greater than 1, and also the lack of xyzw.
'bulk':
    How many of each "bulk type" to scan. A class of constants is
    considered a bulk type iff we expect to have a lot of it in LIReC
    (specified by whether or not it's in the BULK_TYPES array).
'filters':
    A dictionary that specifies which kinds of constants to use to look for relations. Currently supported:
    global filters: 'min_precision' specifies the minimal precision value of constants that will be used.
    'PcfCanonical': 'count' specifies how many pcf's at a time, 'balanced_only' filters to only PCFs of balanced degrees if set to True.
    'Named': 'count' specifies how many named constants at a time.
'''
import mpmath as mp
import os
from functools import reduce
from itertools import combinations, takewhile, product
from logging import getLogger
from logging.config import fileConfig
from sqlalchemy import or_
from sqlalchemy.sql.expression import func
from sympy.polys.polytools import poly_from_expr
from time import time
from traceback import format_exc
from LIReC.config import configuration
from LIReC.lib import models, db_access
from LIReC.lib.pcf import *
from LIReC.lib.pslq_utils import *

mp.mp.dps = 2000

EXECUTE_NEEDS_ARGS = True
DEBUG_PRINT_PRECISION_RATIOS = False
DEBUG_PRINT_CONSTANTS = True

ALGORITHM_NAME = 'POLYNOMIAL_PSLQ'
LOGGER_NAME = 'job_logger'
BULK_SIZE = 500
BULK_TYPES = {'PcfCanonical'}
SUPPORTED_TYPES = ['Named', 'PcfCanonical']
DEFAULT_CONST_COUNT = 1
DEFAULT_DEGREE = 2
DEFAULT_ORDER = 1
MIN_PRECISION_RATIO = 0.8
MAX_PREC = 9999

FILTERS = [
        models.Constant.precision.isnot(None)
        #or_(models.Cf.scanned_algo == None, ~models.Cf.scanned_algo.has_key(ALGORITHM_NAME)) # TODO USE scan_history TABLE!!!
        ]

def get_filters(filters, const_type):
    filter_list = list(FILTERS) # copy!
    if 'global' in filters:
        global_filters = filters['global']
        if 'min_precision' in global_filters:
            filter_list += [models.Constant.precision >= global_filters['min_precision']]
    if const_type == 'PcfCanonical':
        filter_list += [models.PcfCanonicalConstant.convergence != models.PcfConvergence.RATIONAL.value]
        if filters['PcfCanonical'].get('balanced_only', False):
            filter_list += [func.cardinality(models.PcfCanonicalConstant.P) == func.cardinality(models.PcfCanonicalConstant.Q)]

    return filter_list 

def compress_relation(result, consts, exponents, degree, order):
    # will need to use later, so evaluating into lists
    getLogger(LOGGER_NAME).debug(f'Original relation is {result}')
    
    indices_per_var = list(list(i for i, e in enumerate(exponents) if e[j]) for j in range(len(consts)))
    redundant_vars = list(i for i, e in enumerate(indices_per_var) if not any(result[j] for j in e))
    redundant_coeffs = set()
    for redundant_var in sorted(redundant_vars, reverse=True): # remove redundant variables
        getLogger(LOGGER_NAME).debug(f'Removing redundant variable #{redundant_var}')
        redundant_coeffs |= set(indices_per_var[redundant_var])
        del consts[redundant_var]
    
    # remove redundant degrees and orders
    indices_per_degree = list(list(i for i, e in enumerate(exponents) if sum(e.values()) == j) for j in range(degree+1))
    redundant_degrees = list(i for i, e in enumerate(indices_per_degree) if not any(result[j] for j in e))
    redundant_degrees = list(takewhile(lambda x: sum(x) == degree, enumerate(sorted(redundant_degrees, reverse=True))))
    if redundant_degrees:
        degree = redundant_degrees[-1][1] - 1
    redundant_coeffs.update(*indices_per_degree[degree+1:])
    
    indices_per_order = list(list(i for i, e in enumerate(exponents) if max(e.values(), default=0) == j) for j in range(order+1))
    redundant_orders = list(i for i, e in enumerate(indices_per_order) if not any(result[j] for j in e))
    redundant_orders = list(takewhile(lambda x: sum(x) == order, enumerate(sorted(redundant_orders, reverse=True))))
    if redundant_orders:
        order = redundant_orders[-1][1] - 1
    redundant_coeffs.update(*indices_per_order[order+1:])
    
    getLogger(LOGGER_NAME).debug(f'True degree and order are {degree, order}')
    for i in sorted(redundant_coeffs, reverse=True):
        del result[i]
    
    getLogger(LOGGER_NAME).debug(f'Compressed relation is {result}')

    return result, consts, degree, order

def relation_is_new(consts, degree, order, other_relations):
    return not any(r for r in other_relations
                   if {c.const_id for c in r.constants} <= {c.const_id for c in consts} and r.details[0] <= degree and r.details[1] <= order)

def check_consts(consts, exponents, degree, order):
    result, true_prec, min_prec = poly_check(consts, exponents = exponents)
    if not result:
        return []
    if DEBUG_PRINT_PRECISION_RATIOS:
        with mp.workdps(5):
            r = str(true_prec / min_prec)
        print(f'Found relation with precision ratio {r}')
    if true_prec / min_prec < MIN_PRECISION_RATIO:
        if DEBUG_PRINT_PRECISION_RATIOS:
            print(f'Too low! Ignoring...')
        return []
    result, new_consts, new_degree, new_order = compress_relation(result, consts, exponents, degree, order)
    # now must check subrelations! PSLQ is only guaranteed to return a small norm,
    # but not guaranteed to return a 1-dimensional relation, see for example pslq([1,2,3])
    subrelations = []
    for i in range(1, len(consts)):
        exponents = get_exponents(degree, order, i)
        for subset in combinations(consts, i):
            subresult, true_prec2, min_prec2 = poly_check(subset, exponents = exponents)
            if subresult:
                subresult, subconsts, subdegree, suborder = compress_relation(subresult, subset, exponents, degree, order)
                if relation_is_new(subconsts, subdegree, suborder, subrelations) and true_prec2 / min_prec2 >= MIN_PRECISION_RATIO: # no need to check against new_relations + old_relations here btw
                    true_prec2 = min(true_prec2, MAX_PREC)
                    subrelations += [models.Relation(relation_type=ALGORITHM_NAME, details=[subdegree,suborder]+subresult, precision=int(true_prec2), constants=[c.base for c in subconsts])]
    true_prec = min(true_prec, MAX_PREC)
    return subrelations if subrelations else [models.Relation(relation_type=ALGORITHM_NAME, details=[new_degree,new_order]+result, precision=int(true_prec), constants=[c.base for c in new_consts])]

def get_const_class(const_type):
    name = const_type + 'Constant'
    if name not in models.__dict__:
        raise ValueError(f'Unknown constant type {const_type}')
    return models.__dict__[name]

def get_consts_from_query(const_type, query_data):
    const_type = get_const_class(const_type)
    return query_data[[i for i in range(len(query_data)) if isinstance(query_data[i][0], const_type)][0]]

def get_consts(const_type, db, filters):
    if const_type == 'Named': # Constant first intentionally! don't need extra details, but want to filter still
        return db.session.query(models.Constant).join(models.NamedConstant).order_by(models.NamedConstant.const_id).filter(*get_filters(filters, const_type))    

def run_query(filters=None, degree=None, bulk=None):
    fileConfig('LIReC/logging.config', defaults={'log_filename': 'pslq_const_manager'})
    if not filters:
        return []
    bulk_types = set(filters.keys()) & BULK_TYPES
    if not bulk_types:
        return []
    bulk = bulk if bulk else BULK_SIZE
    getLogger(LOGGER_NAME).debug(f'Starting to check relations, using bulk size {bulk}')
    db = db_access.LIReC_DB()
    results = [db.session.query(models.Constant).join(get_const_class(const_type)).filter(*get_filters(filters, const_type)).order_by(func.random()).limit(bulk).all() for const_type in bulk_types]
    # apparently postgresql is really slow with the order_by(random) part,
    # but on 1000 CFs it only takes 1 second, which imo is worth it since
    # that allows us more variety in testing the CFs...
    # TODO what to do if results is unintentionally empty?
    db.session.close()
    getLogger(LOGGER_NAME).info(f'size of batch is {len(results) * bulk}')
    return results

def execute_job(query_data, filters=None, degree=None, order=None, bulk=None, manual=False):
    try: # whole thing must be wrapped so it gets logged
        fileConfig('LIReC/logging.config', defaults={'log_filename': 'analyze_pcfs' if manual else f'pslq_const_worker_{os.getpid()}'})
        global_filters = filters.get('global', {})
        filters.pop('global', 0) # instead of del so we can silently dispose of global even if it doesn't exist
        if not filters:
            getLogger(LOGGER_NAME).error('No filters found! Aborting...')
            return 0 # this shouldn't happen unless pool_handler changes, so just in case...
        keys = filters.keys()
        for const_type in keys:
            if const_type not in SUPPORTED_TYPES:
                msg = f'Unsupported filter type {const_type} will be ignored! Must be one of {SUPPORTED_TYPES}.'
                print(msg)
                getLogger(LOGGER_NAME).warn(msg)
                del filters[const_type]
            elif 'count' not in filters[const_type]:
                filters[const_type]['count'] = DEFAULT_CONST_COUNT
        total_consts = sum(c['count'] for c in filters.values())
        degree = degree if degree else DEFAULT_DEGREE
        order = order if order else DEFAULT_ORDER
        getLogger(LOGGER_NAME).info(f'checking against {total_consts} constants at a time, subdivided into {({k : filters[k]["count"] for k in filters})}, using degree-{degree} relations')
        if degree > total_consts * order:
            degree = total_consts * order
            getLogger(LOGGER_NAME).info(f'redundant degree detected! reducing to {degree}')
        
        db = db_access.LIReC_DB()
        subsets = [combinations(get_consts_from_query(const_type, query_data) if const_type in BULK_TYPES else get_consts(const_type, db, {**filters, 'global':global_filters}), filters[const_type]['count']) for const_type in filters]
        exponents = get_exponents(degree, order, total_consts)
        
        old_relations = db.session.query(models.Relation).all()
        orig_size = len(old_relations)
        # even if the commented code were to be uncommented and implemented for
        # the scan_history table, this loop still can't be turned into list comprehension
        # because finding new relations depends on the new relations we found so far!
        for consts in product(*subsets):
            consts = [c for t in consts for c in t] # need to flatten...
            if relation_is_new(consts, degree, order, old_relations):
                if DEBUG_PRINT_CONSTANTS:
                    getLogger(LOGGER_NAME).debug(f'checking consts: {[c.const_id for c in consts]}')
                new_relations = check_consts(consts, exponents, degree, order)
                if new_relations:
                    getLogger(LOGGER_NAME).info(f'Found relation(s) on constants {[c.const_id for c in consts]}!')
                    try_count = 1
                    while try_count < 3:
                        try:
                            db.session.add_all(new_relations)
                            db.session.commit()
                            old_relations += new_relations
                        except:
                            db.session.rollback()
                            #db.session.close()
                            #db = db_access.LIReC_DB()
                            if try_count == 1:
                                getLogger(LOGGER_NAME).warn('Failed to commit once, trying again.')
                            else:
                                getLogger(LOGGER_NAME).error(f'Could not commit relation(s): {format_exc()}')
                        try_count += 1
            #for cf in consts:
            #    if not cf.scanned_algo:
            #        cf.scanned_algo = dict()
            #    cf.scanned_algo[ALGORITHM_NAME] = int(time())
            #db.session.add_all(consts)
        getLogger(LOGGER_NAME).info(f'finished - found {len(old_relations) - orig_size} results')
        db.session.close()
        
        getLogger(LOGGER_NAME).info('Commit done')
        
        return len(new_relations)
    except:
        getLogger(LOGGER_NAME).error(f'Exception in execute job: {format_exc()}')
        # not returning anything so summarize_results can see the error

def summarize_results(results):
    if not all(results):
        getLogger(LOGGER_NAME).info(f'At least one of the workers had an exception! Check logs')
    getLogger(LOGGER_NAME).info(f'In total found {sum(r for r in results if r)} relations')


PRECISION_LIMIT = 20

def main(): # if testing PCFs that aren't already in LIReC, requires permissions to upload them!
    os.makedirs(os.path.join(os.getcwd(), 'logs'), exist_ok=True)
    config = [c for name, c in configuration if name == 'poly_pslq']
    db = db_access.LIReC_DB()
    pcfs = (PCF(*(poly_from_expr(poly)[0].all_coeffs() for poly in pcf)) for pcf in config['pcfs'])
    successful, unsuccessful = db.add_pcfs(pcfs)
    unsuccessful['Too imprecise'] = []
    for key in unsuccessful:
        if key == 'Already exist': # still want to test these...
            canonical_forms = [[poly.all_coeffs() for poly in pcf.get_canonical_form()] for pcf in unsuccessful[key]]
            for pcf in db.session.query(models.PcfCanonicalConstant):
                if [pcf.P, pcf.Q] in canonical_forms: # TODO theoretically this can be implemented as a filter but i can't figure it out right now...
                    if pcf.base.precision > PRECISION_LIMIT:
                        successful += [pcf]
                    else:
                        unsuccessful['Too imprecise'] += [unsuccessful[key][canonical_forms.index([pcf.P, pcf.Q])]]
        elif unsuccessful[key]:
            polys = [[[int(c) for c in poly.all_coeffs()] for poly in pcf.get_canonical_form()] for pcf in unsuccessful[key]]
            pcfs_string = reduce(lambda a,b: a+'\n\t'+str(b), polys[1:], str(polys[0]))
            print(f'Could not add {len(unsuccessful[key])} pcfs with reason "{key}":\n\t{pcfs_string}\n')
    db.session.close()
    config = config['args']
    if 'PcfCanonical' in config['filters'] and not successful:
        print('None of the configued PCFs can be tested! Aborting...')
    else: # bulk doesn't matter, it just exists due to the limitation of how EXECUTE_NEEDS_ARGS works
        print(f'Testing relations with {len(successful)} PCFs')
        execute_job(successful, config.get('filters', None), config.get('degree', None), config.get('order', None), manual = True)

if __name__ == '__main__':
    '''
    Example config: (using BOINC zetas)
    configuration = {
        ...,
        ('poly_pslq', {
            'args': { ...
            }, ...
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
        })
        ...
    }
    '''
    main()
