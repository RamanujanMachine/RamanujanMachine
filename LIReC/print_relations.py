from functools import reduce
from operator import add
from sympy import Poly, Symbol
import sys
from LIReC.lib.db_access import LIReC_DB
from LIReC.lib.models import *
from LIReC.jobs.job_poly_pslq import get_exponents

def main():
    keep_going = len(sys.argv) > 1
    print(f'printing relations one at a time in descending order of precision{"" if keep_going else ", press enter to print the next"}')
    db = LIReC_DB()
    n = Symbol('n')
    rels = db.session.query(Relation).order_by(Relation.precision.desc()).all()
    nameds = db.session.query(NamedConstant).all()
    pcfs = db.session.query(PcfCanonicalConstant).all()
    for rel in rels:
        exponents = get_exponents(rel.details[:2], len(rel.constants))
        monoms = [reduce(add, (f'*c{i}**{exp[i]}' for i in range(len(rel.constants))), f'{rel.details[2:][j]}') for j, exp in enumerate(exponents)]
        poly = Poly(reduce(add, ['+'+monom for monom in monoms], ''))
        toprint = f'poly: {poly.expr}, precision: {rel.precision}' + ', consts: {\r\n'
        for const in rel.constants:
            named = [n for n in nameds if n.const_id == const.const_id]
            if named:
                toprint += f'    {named[0].name} : {named[0].description}'
            else:
                pcf = [p for p in pcfs if p.const_id == const.const_id]
                if pcf:
                    pcf = pcf[0]
                    if pcf.original_a and pcf.original_b:
                        b_factors = Poly(pcf.original_b, n).factor_list()
                        if len(b_factors[1]) > 1: # usually b_n tends to have only nonpositive integer roots, so try to extract them!
                            b_string = reduce(add, [((f'*({root.expr})' if root.expr != n else '*n') + (f'**{mult}' if mult > 1 else '')) for root, mult in b_factors[1]], f'{b_factors[0]}')
                            toprint += f'    a: {Poly(pcf.original_a, n).expr}, b: {b_string}'
                        else:
                            toprint += f'    a: {Poly(pcf.original_a, n).expr}, b: {Poly(pcf.original_b, n).expr}'
                    else:
                        toprint += f'    P: {Poly(pcf.P, n).expr}, Q: {Poly(pcf.Q, n).expr}'
                else:
                    print(f'constant with uuid {const.const_id} has no known extension')
            toprint += f', precision: {const.precision}, value: {str(const.value)[:50]}...' + '\r\n'
        print(toprint + '}')
        if not keep_going:
            input()

if __name__ == '__main__':
    main()
