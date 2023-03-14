from os import system
from LIReC.config import get_connection_string
from LIReC.lib import models, db_access
from LIReC.lib.calculator import Universal, Constants

if __name__ == '__main__':
    system(f'psql {get_connection_string(db_name="postgres")} < LIReC/lib/create_db.sql')
    
    precision = 4000
    print(f'Using {precision} digits of precision')
    Constants.set_precision(precision)
    db = db_access.LIReC_DB()
    for const in Constants.__dict__.keys():
        if const[0] == '_' or const == 'set_precision':
            continue
        print(f'Adding named constant {const}')
        db.session.add(Universal.calc_named(const, True, True))
    
    from mpmath import zeta
    
    for x in [2, 4, 5, 6, 7]:
        named_const = models.NamedConstant()
        named_const.base = models.Constant()
        named_const.base.precision = precision
        named_const.base.value = Decimal(str(zeta(x)))
        named_const.name = f'Zeta{x}'
        named_const.description = f'zeta({x})'
        db.session.add(named_const)
    
    db.session.commit()
    db.session.close()
