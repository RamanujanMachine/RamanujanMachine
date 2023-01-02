from LIReC.lib.db_access import LIReC_DB
from LIReC.lib.models import *

_db = LIReC_DB()

# >>> Remember to first configure config.py with your username+password! <<<

def next() -> PcfCanonicalConstant or None:
    '''
    Returns the highest priority oldest untweeted PCF.
    '''
    try:
        return _db.session.query(PcfCanonicalConstant).join(Constant)
                  .order_by(Constant.priority.desc(), Constant.time_added)
                  .filter(Constant.tweeted == 0).first()
    except Exception as e:
        _db.session.rollback()
        raise e

def next_relation(const_id: str) -> Relation or None:
    '''
    Returns the highest priority oldest untweeted relation in which the given constant participates.
    '''
    try:
        return _db.session.query(Relation)
                          .order_by(Constant.priority.desc(), Constant.time_added)
                          .filter(Relation.constants.any(Constant.const_id == const_id)).first()
    except Exception as e:
        _db.session.rollback()
        raise e

def set_tweeted(obj: Constant or Relation, auto_commit : bool = True) -> None:
    '''
    Updates the given object to indicate it has been tweeted. This will cause it to be ignored in future calls to next() or next_relation(const_id).
    '''
    obj.tweeted = 1
    if auto_commit:
        try:
            _db.session.commit()
        except Exception as e:
            _db.session.rollback()
            raise e
