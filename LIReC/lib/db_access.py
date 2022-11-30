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
from LIReC.config import get_connection_string, configuration
from LIReC.lib import models
from LIReC.lib.pcf import *

class LIReC_DB:
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
                                calculation : PCFCalc or None = None) -> models.PcfCanonicalConstant:
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

    def add_pcf(self, pcf: PCF) -> None:
        """
        Expect PCF object.
        raises IntegrityError if pcf already exists in LIReC.
        raises NoFRException if calculate is True and the pcf doesn't converge
        raises IllegalPCFException if the pcf has natural roots or if its b_n has irrational roots.
        """
        if any(r for r in pcf.a.real_roots() if isinstance(r, Integer) and r > 0):
            raise PcfCalc.IllegalPCFException('Natural root in partial denominator ensures divergence.')
        if any(r for r in pcf.b.real_roots() if isinstance(r, Integer) and r > 0):
            raise PcfCalc.IllegalPCFException('Natural root in partial numerator ensures trivial convergence to a rational number.')
        if any(r for r in pcf.b.all_roots() if not isinstance(r, Rational)):
            raise PcfCalc.IllegalPCFException('Irrational or Complex roots in partial numerator are not allowed.')
        top, bot = pcf.get_canonical_form()
        #calculation = LIReC_DB.calc_pcf(pcf, depth) if depth else None
        # By default the coefs are sympy.core.numbers.Integer but sql need them to be integers
        return self.add_pcf_canonical([int(coef) for coef in top.all_coeffs()], [int(coef) for coef in bot.all_coeffs()], PcfCalc(pcf).run(**configuration['auto_pcf']))
    
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
                    raise e # otherwise already in LIReC
                unsuccessful['Already exist'] += [pcf]
            except PcfCalc.NoFRException:
                unsuccessful['No FR'] += [pcf]
            except PcfCalc.IllegalPCFException:
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
                    raise e # otherwise already in LIReC
            except PcfCalc.NoFRException:
                pass
            except PcfCalc.IllegalPCFException:
                pass
    
    @staticmethod
    def parse_cf_to_lists(cf: models.PcfCanonicalConstant) -> CanonicalForm:
        return [int(coef) for coef in cf.P], [int(coef) for coef in cf.Q]
    
    def get_canonical_forms(self) -> List[CanonicalForm]:
        return [LIReC_DB.parse_cf_to_lists(pcf) for pcf in self.cfs.all()]
    
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
