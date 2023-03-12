'''
Gradually "refines" constants in LIReC by taking the constants with lowest precision,
increasing their precision, and then rechecking the relations they are participating in.

Configured as such:
    'consts_per_core': Number of constants to refine per core (per iteration). Defaults to 5.
    'greedy_precision': If nonzero, will instead refine constants whose precision is just below greedy_precision (if there are enough). Defaults to 0.
'''
from logging import getLogger
from logging.config import fileConfig
from os import getpid
from functools import reduce
from operator import add
from os import cpu_count
from traceback import format_exc
from LIReC.config import configuration
from LIReC.lib import models, db_access
from LIReC.lib.calculator import Universal

LOGGER_NAME = 'job_logger'
EXTENSION_TYPES = [models.NamedConstant, models.DerivedConstant, models.PcfCanonicalConstant] # ordered by rough preference to compute (pcfs can be slow)
DEBUG_PRINT = True

def run_query(consts_per_core=5, greedy_precision=0):
    fileConfig('LIReC/logging.config', defaults={'log_filename': 'refining_manager'})
    total_consts = configuration['jobs_to_run']['refine'].get('async_cores', cpu_count()) * consts_per_core
    db = db_access.LIReC_DB()
    results = []
    while len(results) < total_consts:
        query = db.session.query(models.Constant)
        if greedy_precision:
            query = query.filter(models.Constant.precision < greedy_precision).order_by(models.Constant.precision.desc())
        else:
            query = query.order_by(models.Constant.precision)
        results = query.limit(total_consts).all()
        if greedy_precision and len(results) < total_consts: # not enough!
            greedy_precision = 0
    db.session.close()
    getLogger(LOGGER_NAME).info(f'size of batch is {len(results)}')
    if DEBUG_PRINT:
        getLogger(LOGGER_NAME).debug('original values are:')
        for const in results:
            getLogger(LOGGER_NAME).debug(f'{const.const_id}: precision {const.precision}, value {const.value}')
    return results

def get_extensions(db, consts):
    ids = [c.const_id for c in consts]
    by_type = [db.session.query(ext).filter(ext.const_id.in_(ids)).all() for ext in EXTENSION_TYPES]
    return [reduce(add, [[c for c in ext if c.const_id == const_id] for ext in by_type]) for const_id in ids]

def execute_job(query_data):
    try: # whole thing must be wrapped so it gets logged
        fileConfig('LIReC/logging.config', defaults={'log_filename': f'refining_worker_{getpid()}'})
        db = db_access.LIReC_DB()
        extensions = get_extensions(db, query_data)
        refined = []
        for const, exts in zip(query_data, extensions):
            if DEBUG_PRINT:
                getLogger(LOGGER_NAME).debug(f'refining constant #{const.const_id}')
            Universal.set_precision(const.precision * 2) # will be ignored for every ext in exts that doesn't rely on precision, like PcfCanonicalConstant which uses depth instead
            res = None
            for ext in exts:
                res = Universal.calc_silent(ext, db, const)
                if res:
                    refined += [res]
                    break
            if not res: # communicates to central job that no refining happened
                refined += [None]
        
        if DEBUG_PRINT:
            getLogger(LOGGER_NAME).debug(f'done')
        return refined
    except:
        getLogger(LOGGER_NAME).error(f'Exception in execute job: {format_exc()}')

def summarize_results(results):
    db = db_access.LIReC_DB()
    try:
        results = [r for r in results if r]
        getLogger(LOGGER_NAME).info(f'refined a total of {len(results)} constants')
        if DEBUG_PRINT:
            getLogger(LOGGER_NAME).debug('new values are:')
            for const in results:
                getLogger(LOGGER_NAME).debug(f'{const.const_id}: precision {const.base.precision}, value {const.base.value}')
        db.session.add_all(results)
        db.session.commit()
    except:
        getLogger(LOGGER_NAME).error(f'Exception while trying to commit refines: {format_exc()}')
        db.session.rollback()
    
    db.session.close()
