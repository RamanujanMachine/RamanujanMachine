import re
from decimal import Decimal
import os
from traceback import format_exc
from urllib.request import urlopen
from LIReC.config import get_connection_string
from LIReC.lib import models, db_access
from LIReC.lib.calculator import Constants

COMMAND = 'psql {connection_string} < LIReC/lib/create_db.sql'
MIN_PRECISION = 20
PREFIX_LINK = 'OEIS link: '
PREFIX_URL = 'org/A'
TABLE_FORMAT = '{0}/b{1}.txt'
FORMATTING = 'utf-8'

def read_oeis(file): # works with the OEIS format, credit where credit is due
    value = ''
    first_index = None
    precision = None
    while True:
        line = file.readline().decode(FORMATTING)
        if not line:
            return value, int(precision)
        res = re.match('\\s*(\\d+)\\s*(\\d+)\\s*', line)
        if not res:
            continue
        if not value:
            first_index = int(res.group(1)) # this is the number of digits before the decimal point!
        if first_index == 0:
            if not value:
                value = '0'
            value += '.'
        value += res.group(2)
        first_index -= 1
        precision = res.group(1)

if __name__ == '__main__':
    os.system(COMMAND.format(connection_string=get_connection_string(db_name='postgres')))
    
    precision = 4000
    print(f'Using {precision} digits of precision')
    Constants.set_precision(precision)
    db = db_access.LIReC_DB()
    for const in Constants.__dict__.keys():
        if const[0] == '_' or const == 'set_precision':
            continue
        print(f'Adding named constant {const}')
        named_const = models.NamedConstant()
        const_func = Constants.__dict__[const].__get__(0) # no idea why the 0 here is needed, but it's needed alright...
        named_const.base = models.Constant()
        if 'WARNING' not in const_func.__doc__: # too slow to calculate!!
            if 'CAUTION' in const_func.__doc__:
                print(f'    Calculation of {const} is expected to take somewhat longer...')
            named_const.base.precision = precision
            named_const.base.value = Decimal(str(const_func()))
        else:
            print(f'    Skipping calculation of {const}, too inefficient or no calculation available!')
            if PREFIX_LINK not in const_func.__doc__:
                print(f'No backup value exists for {const}! Add one when possible (probably using OEIS)')
            else:
                i = const_func.__doc__.index(PREFIX_LINK)
                url = const_func.__doc__[i + len(PREFIX_LINK) : i + const_func.__doc__[i:].index('\n')]
                url = TABLE_FORMAT.format(url, url[url.index(PREFIX_URL) + len(PREFIX_URL) : ])
                try:
                    value, precision2 = read_oeis(urlopen(url))
                    if precision2 < MIN_PRECISION:
                        print(f'    OEIS value has too low precision, check back if and when it has at least {MIN_PRECISION} digits.')
                    else:
                        print('    OEIS value found, will be used instead')
                        named_const.base.value = value[:16001] # the numeric type is limited to 16383 digits after the decimal point apparently, so for now this sits here
                        named_const.base.precision = min(precision2, 16000)
                except:
                    print(f'Exception while fetching {const} from OEIS: {format_exc()}')
                
        named_const.name = const
        named_const.description = const_func.__doc__[:const_func.__doc__.index('.\n')]
        db.session.add(named_const)
    
    from mpmath import zeta
    
    for x in [2, 4, 5, 6, 7]:
        named_const = models.NamedConstant()
        named_const.base = models.Constant()
        named_const.base.precision = 4000
        named_const.base.value = Decimal(str(zeta(x)))
        named_const.name = f'Zeta{x}'
        named_const.description = f'zeta({x})'
        db.session.add(named_const)
    
    db.session.commit()
    db.session.close()
