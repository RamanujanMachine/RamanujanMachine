from logging import getLogger
from logging.config import fileConfig
from traceback import format_exc
from LIReC.lib import models, db_access
from LIReC.lib.pslq_utils import poly_verify

LOGGER_NAME = 'job_logger'
PRECISION_TOLERANCE = 1.1

def execute_job():
    try:
        fileConfig('logging.config', defaults={'log_filename': 'verify'})
        db = db_access.LIReC_DB()
        rels = db.session.query(models.Relation).all()
        
        # way faster to query in bulk and filter locally than to query relations individually!
        consts = db.session.query(models.Constant).all()
        consts = {c.const_id : c for c in consts}
        mtm = db.session.query(models.constant_in_relation_table).all()
        rels = [(r, [consts[row[0]] for row in mtm if row[1] == r.relation_id]) for r in rels]
        for r, constants in rels:
            true_precision = poly_verify(constants, full_relation = r.details)
            if true_precision * PRECISION_TOLERANCE >= r.precision:
                r.precision = true_precision
            else: # too bad! don't over-report your results next time
                db.session.delete(r)
        db.session.commit()
        db.session.close()
    except:
        getLogger(LOGGER_NAME).error(f'Exception while verifying: {format_exc()}')

def main():
    execute_job()

if __name__ == '__main__':
    main()
