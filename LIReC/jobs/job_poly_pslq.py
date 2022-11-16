'''
Finds polynomial relations between constants in LIReC, using PSLQ.

Configured as such:
'degree': Tuple of the form (polydegree: int, innerdegree: int). All relations are structured like
          multivariate polynomials over the constants and CFs, of degree polydegree with a maximum
          exponent of innerdegree. For example, a 2-variable polynomial of degree (2,1) will be of
          the form a+bx+cy+dxy (note the lack of x^2 and y^2), and a 4-variable polynomial of degree
          (3,1) will be of the form:
              a + bx+cy+dz+ew + fxy+gxz+hxw+iyz+jyw+kzw + lxyz+mxyw+nxzw+oyzw
          Note here the lack of any single variable with an exponent greater than 1, and also the lack of xyzw.
'bulk': How many of each "bulk type" to scan. A class of constants is
        considered a bulk type iff we expect to have a lot of it in LIReC
        (specified by whether or not it's in the BULK_TYPES array).
'subdivide': A dictionary that specifies which kinds of constants to use to look for relations. Currently supported:
    'PcfCanonical': 'count' specifies how many pcf's at a time, 'balanced_only' filters to only PCFs of balanced degrees if set to True.
    'Named': 'count' specifies how many named constants at a time.
'''
import mpmath as mp
from time import time
from sqlalchemy import Integer, or_, Float
from sqlalchemy.orm.attributes import flag_modified
from sqlalchemy.sql.expression import func
from logging import getLogger
from logging.config import fileConfig
from os import getpid
from itertools import chain, combinations, combinations_with_replacement, takewhile, product
from collections import Counter
from functools import reduce
from operator import mul
from traceback import format_exc
from LIReC.lib import models, db_access

mp.mp.dps = 2000

EXECUTE_NEEDS_ARGS = True
DEBUG_PRINT_PRECISION_RATIOS = True

ALGORITHM_NAME = 'POLYNOMIAL_PSLQ'
LOGGER_NAME = 'job_logger'
BULK_SIZE = 500
BULK_TYPES = {'PcfCanonical'}
SUPPORTED_TYPES = ['Named', 'PcfCanonical']
DEFAULT_CONST_COUNT = 1
DEFAULT_DEGREE = (2, 1)
MIN_PRECISION_RATIO = 0.8
MAX_PREC = 9999

FILTERS = [
        models.Constant.precision.isnot(None)
        #or_(models.Cf.scanned_algo == None, ~models.Cf.scanned_algo.has_key(ALGORITHM_NAME)) # TODO USE scan_history TABLE!!!
        ]

def get_filters(subdivide, const_type):
    filters = list(FILTERS) # copy!
    if const_type == 'PcfCanonical':
        filters += [models.PcfCanonicalConstant.convergence != models.PcfConvergence.RATIONAL.value]
        if subdivide['PcfCanonical'].get('balanced_only', False):
            filters += [func.cardinality(models.PcfCanonicalConstant.P) == func.cardinality(models.PcfCanonicalConstant.Q)]

    return filters 

def poly_check(consts, exponents):
    mp.mp.dps = min(c.base.precision for c in consts)
    values = [mp.mpf(str(c.base.value)) for c in consts]
    if 1 in values: # solely for backwards-compatibility. We don't need 1 in the DB!
        return None, None, None
    poly = [reduce(mul, (values[i] ** exp[i] for i in range(len(values))), mp.mpf(1)) for exp in exponents]
    try:
        mp.mp.dps = 15 # intentionally low-resolution to quickly try something basic...
        res = mp.pslq(poly)
        if res: # then calculating the substance in what we just found!
            mp.mp.dps = max(c.base.precision for c in consts) + 10
            min_prec = min(c.base.precision for c in consts)
            return res, min(mp.floor(-mp.log10(abs(mp.fdot(poly, res)))), min_prec), min_prec
    except ValueError:
        # one of the constants has too small precision, or one constant
        # is small enough that another constant is smaller than its precision.
        # eitherway there's no relation to be found here!
        pass 
    return None, None, None

def compress_relation(result, consts, exponents, degree):
    # will need to use later, so evaluating into lists
    getLogger(LOGGER_NAME).info(f'Original relation is {result}')
    
    indices_per_var = list(list(i[0] for i in enumerate(exponents) if i[1][j]) for j in range(len(consts)))
    redundant_vars = list(i[0] for i in enumerate(indices_per_var) if not any(result[j] for j in i[1]))
    redundant_coeffs = set()
    for redundant_var in redundant_vars: # remove redundant variables
        getLogger(LOGGER_NAME).info(f'Removing redundant variable #{redundant_var}')
        redundant_coeffs |= set(indices_per_var[redundant_var])
        consts = consts[:redundant_var] + consts[redundant_var + 1:]
    
    polydegree, innerdegree = degree # remove redundant degrees
    indices_per_polydegree = list(list(i[0] for i in enumerate(exponents) if sum(i[1].values())==j) for j in range(polydegree+1))
    redundant_polydegrees = list(i[0] for i in enumerate(indices_per_polydegree) if not any(result[j] for j in i[1]))
    redundant_polydegrees = list(takewhile(lambda x: sum(x) == polydegree, enumerate(sorted(redundant_polydegrees, reverse=True))))
    if redundant_polydegrees:
        polydegree = redundant_polydegrees[-1][1] - 1
    redundant_coeffs.update(*indices_per_polydegree[polydegree+1:])
    
    indices_per_innerdegree = list(list(i[0] for i in enumerate(exponents) if max(i[1].values(), default=0)==j) for j in range(innerdegree+1))
    redundant_innerdegrees = list(i[0] for i in enumerate(indices_per_innerdegree) if not any(result[j] for j in i[1]))
    redundant_innerdegrees = list(takewhile(lambda x: sum(x) == innerdegree, enumerate(sorted(redundant_innerdegrees, reverse=True))))
    if redundant_innerdegrees:
        innerdegree = redundant_innerdegrees[-1][1] - 1
    redundant_coeffs.update(*indices_per_innerdegree[innerdegree+1:])
    
    degree = [polydegree, innerdegree]
    getLogger(LOGGER_NAME).info(f'True degree is {degree}')
    for i in sorted(redundant_coeffs, reverse=True):
        del result[i]
    
    getLogger(LOGGER_NAME).info(f'Compressed relation is {result}')

    return result, consts, degree

def get_exponents(degree, total_consts):
    polydegree, innerdegree = degree
    return list(c for c in map(Counter, chain.from_iterable(combinations_with_replacement(range(total_consts), i) for i in range(polydegree+1)))
                if not any(i for i in c.values() if i > innerdegree))

def relation_is_new(consts, degree, other_relations):
    return not any(r for r in other_relations
                   if {c.const_id for c in r.constants} <= {c.const_id for c in consts} and r.details[0] <= degree[0] and r.details[1] <= degree[1])

def check_consts(consts, exponents, degree):
    result, true_prec, min_prec = poly_check(consts, exponents)
    if not result:
        return []
    with mp.workdps(5):
        r = str(true_prec / min_prec)
    if DEBUG_PRINT_PRECISION_RATIOS:
        print(f'Found relation with precision ratio {r}')
    if true_prec / min_prec < MIN_PRECISION_RATIO:
        if DEBUG_PRINT_PRECISION_RATIOS:
            print(f'Too low! Ignoring...')
        return []
    result, new_consts, new_degree = compress_relation(result, consts, exponents, degree)
    # now must check subrelations! PSLQ is only guaranteed to return a small norm,
    # but not guaranteed to return a 1-dimensional relation, see for example pslq([1,2,3])
    subrelations = []
    for i in range(1, len(consts)):
        exponents = get_exponents(degree, i)
        for subset in combinations(consts, i):
            subresult, true_prec2, min_prec2 = poly_check(subset, exponents)
            if subresult:
                subresult, subconsts, subdegree = compress_relation(subresult, subset, exponents, degree)
                if relation_is_new(subconsts, subdegree, subrelations) and true_prec2 / min_prec2 >= MIN_PRECISION_RATIO: # no need to check against new_relations + old_relations here btw
                    true_prec2 = min(true_prec2, MAX_PREC)
                    subrelations += [models.Relation(relation_type=ALGORITHM_NAME, details=subdegree+subresult, precision=int(true_prec2), constants=[c.base for c in subconsts])]
    true_prec = min(true_prec, MAX_PREC)
    return subrelations if subrelations else [models.Relation(relation_type=ALGORITHM_NAME, details=new_degree+result, precision=int(true_prec), constants=[c.base for c in new_consts])]

def transpose(l):
    return [[l[j][i] for j in range(len(l))] for i in range(len(l[0]))] if l else [[]]

def get_consts_from_query(const_type, query_data):
    return query_data[[i for i in range(len(query_data)) if isinstance(query_data[i][0], eval(f'models.{const_type}Constant'))][0]]

def get_consts(const_type, db, subdivide):
    if const_type == 'Named':
        return db.constants.join(models.Constant).filter(models.Constant.value.isnot(None))

def run_query(subdivide=None, degree=None, bulk=None):
    fileConfig('LIReC/logging.config', defaults={'log_filename': 'pslq_const_manager'})
    if not subdivide:
        return []
    bulk_types = set(subdivide.keys()) & BULK_TYPES
    if not bulk_types:
        return []
    bulk = bulk if bulk else BULK_SIZE
    getLogger(LOGGER_NAME).debug(f'Starting to check relations, using bulk size {bulk}')
    db= db_access.LIReC_DB()
    results = [db.session.query(eval(f'models.{const_type}Constant')).join(models.Constant).filter(*get_filters(subdivide, const_type)).order_by(func.random()).limit(bulk).all() for const_type in bulk_types]
    # apparently postgresql is really slow with the order_by(random) part,
    # but on 1000 CFs it only takes 1 second, which imo is worth it since
    # that allows us more variety in testing the CFs...
    # TODO what to do if results is unintentionally empty?
    db.session.close()
    getLogger(LOGGER_NAME).info(f'size of batch is {len(results) * bulk}')
    return transpose(results) # so pool_handler can correctly divide among the sub-processes

def execute_job(query_data, subdivide=None, degree=None, bulk=None, manual=False):
    fileConfig('LIReC/logging.config', defaults={'log_filename': 'analyze_pcfs' if manual else f'pslq_const_worker_{getpid()}'})
    if not subdivide:
        getLogger(LOGGER_NAME).error('Nothing to do! Aborting...')
        return 0 # this shouldn't happen unless pool_handler changes, so just in case...
    keys = subdivide.keys()
    for const_type in keys:
        if const_type not in SUPPORTED_TYPES:
            msg = f'Unsupported constant type {const_type} will be ignored! Must be one of {SUPPORTED_TYPES}.'
            print(msg)
            getLogger(LOGGER_NAME).warn(msg)
            del subdivide[const_type]
        elif 'count' not in subdivide[const_type]:
            subdivide[const_type]['count'] = DEFAULT_CONST_COUNT
    total_consts = sum(c['count'] for c in subdivide.values())
    degree = degree if degree else DEFAULT_DEGREE
    getLogger(LOGGER_NAME).info(f'checking against {total_consts} constants at a time, subdivided into {({k : subdivide[k]["count"] for k in subdivide})}, using degree-{degree} relations')
    if degree[0] > total_consts * degree[1]:
        degree = (total_consts * degree[1], degree[1])
        getLogger(LOGGER_NAME).info(f'redundant degree detected! reducing to {degree}')
    
    query_data = transpose(query_data)
    try:
        db = db_access.LIReC_DB()
        subsets = [combinations(get_consts_from_query(const_type, query_data) if const_type in BULK_TYPES else get_consts(const_type, db, subdivide), subdivide[const_type]['count']) for const_type in subdivide]
        exponents = get_exponents(degree, total_consts)
        
        old_relations = db.session.query(models.Relation).all()
        orig_size = len(old_relations)
        # even if the commented code were to be uncommented and implemented for
        # the scan_history table, this loop still can't be turned into list comprehension
        # because finding new relations depends on the new relations we found so far!
        for consts in product(*subsets):
            consts = [c for t in consts for c in t] # need to flatten...
            if relation_is_new(consts, degree, old_relations):
                getLogger(LOGGER_NAME).info(f'checking consts: {[c.const_id for c in consts]}')
                new_relations = check_consts(consts, exponents, degree)
                if new_relations:
                    getLogger(LOGGER_NAME).info(f'Found relation(s)!')
                    old_relations += new_relations
                    db.session.add_all(new_relations)
                    db.session.commit()
                    
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

def summarize_results(results):
    getLogger(LOGGER_NAME).info(f'In total found {sum(results)} relations')
