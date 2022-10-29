from __future__ import annotations
import logging
import mpmath as mp
import numpy as np
from collections import namedtuple
from decimal import Decimal, getcontext
from functools import reduce
from sympy import Poly
from sympy.core.numbers import Integer, Rational
from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError
from sqlalchemy.dialects import postgresql
from sqlalchemy.orm import sessionmaker
from psycopg2.errors import UniqueViolation
from typing import Tuple, List, Dict, Generator
from db.config import get_connection_string
from db.lib import models
from db.lib.pcf import *

MAX_ITERATOR_STEPS = 1402
DEPTH = 12000
CALC_JUMP = 200
FR_THRESHOLD = 0.1

PcfCalculation = namedtuple('PcfCalculation', ['value', 'precision', 'last_matrix', 'depth', 'convergence'])

getcontext().prec = 2000
mp.mp.dps = 2000

class PcfCalc:
    a: Poly
    b: Poly
    mat: mp.matrix
    depth: int
    
    def __init__(self: PcfCalc, a: Poly, b: Poly, prev: List[int] or None = None, depth: int = 0):
        self.a = a
        self.b = b
        #self.reduction = 1 # used to exist in the old code, not sure what use we have for this
        self.mat = mp.matrix([prev[0:2], prev[2:4]] if prev else [[a(0), 1], [1, 0]])
        self.depth = depth
    
    def refine(self: PcfCalc):
        # TODO upgrade to mergesort-like binary tree scheme?
        self.depth += 1
        mat = mp.matrix([[self.a(self.depth), self.b(self.depth)], [1, 0]]) * self.mat
        # the old code only did this gcd step every REDUCE_JUMP steps, for now I
        # decided to do it every step to see how much of a performance hit it really is
        gcd = np.gcd.reduce(self.mat)
        #self.reduction *= gcd
        self.mat = mat / mp.mpf(gcd)
    
    @property
    def value(self: PcfCalc):
        return self.mat[0,0] / self.mat[0,1]
    
    @property
    def precision(self: PcfCalc):
        return mp.floor(-mp.log10(abs(self.value - self.mat[1,0] / self.mat[1,1]))) if all(self.mat[:,1]) else -mp.inf

class RamanujanDB(object):
    def __init__(self):
        logging.debug("Trying to connect to database")
        self._engine = create_engine(get_connection_string(), echo=False)
        Session = sessionmaker(bind=self._engine)
        self.session = Session() # TODO can be unit-tested with mock-alchemy or alchemy-mock...
        logging.debug("Connected to database")

    @property
    def constants(self):
        return self.session.query(models.NamedConstant).order_by(models.NamedConstant.const_id)

    @property
    def cfs(self):
        return self.session.query(models.PcfCanonicalConstant).order_by(models.PcfCanonicalConstant.const_id)

    def add_pcf_canonical(self, canonical_form_numerator: List[int], canonical_form_denominator: List[int],
                                calculation : PcfCalculation or None = None) -> models.PcfCanonicalConstant:
        # TODO implement add_pcf_canonicals that uploads multiple at a time
        pcf = models.PcfCanonicalConstant()
        pcf.base = models.Constant()
        pcf.P = canonical_form_numerator
        pcf.Q = canonical_form_denominator
        if calculation:
            pcf.base.value = calculation.value
            pcf.base.precision = calculation.precision
            pcf.last_matrix = reduce(lambda a,b: a+','+str(b), calculation.last_matrix[1:], str(calculation.last_matrix[0]))
            pcf.depth = calculation.depth
            pcf.convergence = calculation.convergence
        self.session.add(pcf)
        
        # yes, commit and check error is better than preemptively checking if unique and then adding,
        # since the latter is two SQL commands instead of one, which breaks on "multithreading" for example
        # and also should be generally slower
        # also this can't be turned into a kind of INSERT pcf ON CONFLICT DO NOTHING statement
        # since this needs the base Constant to be added first so it gains its const_id...
        # TODO investigate if triggers can make something like ON CONFLICT DO NOTHING work anyway,
        # possibly will help with the previous TODO... maybe something like:  https://stackoverflow.com/questions/46105982/postgres-trigger-function-on-conflict-update-another-table
        try:
            self.session.commit()
            return pcf
        except Exception as e:
            self.session.rollback()
            raise e
    
    class IllegalPCFException(Exception):
        pass

    class NoFRException(Exception):
        pass
    
    @staticmethod
    def check_convergence(p, q, fr_list) -> models.PcfConvergence:
        val = p / q
        if mp.almosteq(p, 0) or mp.almosteq(val, 0):
            logging.debug(f'Checking rational that is too close to 0 p={p},q={q}')
            return models.PcfConvergence.RATIONAL
        
        if mp.pslq([val, 1], tol=mp.power(10, -100)):
            return models.PcfConvergence.RATIONAL
        
        if any(abs(fr_list[i + 1] - fr_list[i]) < FR_THRESHOLD for i in range(len(fr_list) - 1)):
            return models.PcfConvergence.FR
        
        if any(abs(fr_list[i + 1] - fr_list[i + 2]) > abs(fr_list[i] - fr_list[i + 1]) for i in range(len(fr_list) - 2)):
            raise RamanujanDB.NoFRException()
        
        return models.PcfConvergence.INDETERMINATE_FR
    
    @staticmethod
    def calc_pcf(pcf):
        # P.S.: The code here is in fact similar to enumerators.FREnumerator.check_for_fr, but
        # the analysis we're doing here is both more delicate (as we allow numerically-indeterminate PCFs),
        # and also less redundant (since we also want the value of the PCF instead of just discarding it for instance)
        calc = PcfCalc(pcf.a, pcf.b)
        
        fr_list = []
        for n in range(1, DEPTH):
            calc.refine()
            if n % CALC_JUMP == 0:
                # yes this calculation converges as i -> inf if there's factorial
                # reduction, ask me for proof if you can't figure out the details. (Itay)
                # (though as far as i could tell, the reverse direction might not be necessarily true)
                fr_list += [mp.log(mp.mpf(np.gcd(*calc.mat[0,:]))) / mp.mpf(n) + calc.a.degree() * (1 - mp.log(n))]
        
        prec = calc.precision
        if prec == -mp.inf:
            raise RamanujanDB.IllegalPCFException('continuant denominator zero')
        
        value = calc.value
        if value and mp.almosteq(0, value):
            logging.debug('Rounding to 0')
            value = 0
        return PcfCalculation(Decimal(str(value)), 2000 if prec == mp.inf else int(prec), [int(x) for x in calc.mat], DEPTH, RamanujanDB.check_convergence(*calc.mat[0,:], fr_list).value)

    def add_pcf(self, pcf: PCF, calculate: bool = True) -> None:
        """
        Expect PCF object.
        raises IntegrityError if pcf already exists in the db.
        raises NoFRException if calculate is True and the pcf doesn't converge
        raises IllegalPCFException if the pcf has natural roots or if its b_n has irrational roots.
        """
        if any(r for r in pcf.a.real_roots() if isinstance(r, Integer) and r > 0):
            raise RamanujanDB.IllegalPCFException('Natural root in partial denominator ensures divergence.')
        if any(r for r in pcf.b.real_roots() if isinstance(r, Integer) and r > 0):
            raise RamanujanDB.IllegalPCFException('Natural root in partial numerator ensures trivial convergence to a rational number.')
        if any(r for r in pcf.b.all_roots() if not isinstance(r, Rational)):
            raise RamanujanDB.IllegalPCFException('Irrational or Complex roots in partial numerator are not allowed.')
        top, bot = pcf.get_canonical_form()
        calculation = RamanujanDB.calc_pcf(pcf) if calculate else None
        # By default the coefs are sympy.core.numbers.Integer but sql need them to be integers
        return self.add_pcf_canonical([int(coef) for coef in top.all_coeffs()], [int(coef) for coef in bot.all_coeffs()], calculation)
    
    def add_pcfs(self, pcfs: Generator[PCF, None, None]) -> Tuple[List[models.PcfCanonicalConstant], Dict[str, List[PCF]]]:
        """
        Expects a list of PCF objects.
        """
        successful = []
        unsuccessful = {'Already exist': [], 'No FR': [], 'Illegal': []}
        for pcf in pcfs:
            try:
                successful.append(self.add_pcf(pcf))
            except IntegrityError as e:
                if not isinstance(e.orig, UniqueViolation):
                    raise e # otherwise already in the DB
                unsuccessful['Already exist'] += [pcf]
            except RamanujanDB.NoFRException:
                unsuccessful['No FR'] += [pcf]
            except RamanujanDB.IllegalPCFException:
                unsuccessful['Illegal'] += [pcf]
        return successful, unsuccessful
    
    def add_pcfs_silent(self, pcfs: Generator[PCF, None, None]) -> None:
        """
        Expects a list of PCF objects. Doesn't return which PCFs were successfully or unsuccessfully added.
        """
        for pcf in pcfs:
            try:
                self.add_pcf(pcf)
            except IntegrityError as e:
                if not isinstance(e.orig, UniqueViolation):
                    raise e # otherwise already in the DB
            except RamanujanDB.NoFRException:
                pass
            except RamanujanDB.IllegalPCFException:
                pass
    
    @staticmethod
    def parse_cf_to_lists(cf: models.PcfCanonicalConstant) -> CanonicalForm:
        return [int(coef) for coef in cf.P], [int(coef) for coef in cf.Q]
    
    def get_canonical_forms(self) -> List[CanonicalForm]:
        return [RamanujanDB.parse_cf_to_lists(pcf) for pcf in self.cfs.all()]
    
    def get_actual_pcfs(self) -> List[PCF]:
        """
        return a list of PCFs from the DB
        """
        return [PCF.from_canonical_form(c) for c in self.get_canonical_forms()]


def main():
    [print(pcf) for pcf in get_actual_pcfs_from_db()]
    print("aa")


if __name__ == "__main__":
    main()
