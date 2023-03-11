from decimal import Decimal, getcontext
from functools import reduce
import numpy as np
import mpmath as mp
from mpmath import mpf
from operator import mul
from re import match
from sympy import Poly, Symbol, compose
from typing import List
from urllib.request import urlopen
from LIReC.lib.models import *
from LIReC.lib.pcf import *

class Universal:
    '''
    central hub for calculating any constant from lib.models to arbitrary precision (if possible).
    '''
    
    @staticmethod
    def set_precision(prec: int = 4000) -> None:
        '''
        set the precision (in significant digits in base 10). shared with other classes in this file
        '''
        mp.mp.dps = prec
    
    @staticmethod
    def read_oeis(file): # works with the OEIS format, credit where credit is due
        value = ''
        first_index = None
        precision = None
        while True:
            line = file.readline().decode('utf-8')
            if not line:
                return value, int(precision)
            res = match(r'\s*(\d+)\s*(\d+)\s*', line)
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
    
    @staticmethod
    def calc_named(name: str or NamedConstant, base=None, verbose=False, force=False):
        PREFIX_LINK = 'OEIS link: '
        PREFIX_URL = 'org/A'
        MIN_PRECISION = 20
        
        const = name if isinstance(name, str) else name.name
        if const not in Constants.__dict__:
            return None # who dis
        named_const = name if isinstance(name, NamedConstant) else NamedConstant()
        const_func = Constants.__dict__[const].__get__(0) # no idea why the 0 here is needed, but it's needed alright...
        named_const.base = base if base else Constant()
        if 'WARNING' not in const_func.__doc__: # too slow to calculate!!
            if verbose and 'CAUTION' in const_func.__doc__:
                print(f'    Calculation of {const} is expected to take somewhat longer...')
            named_const.base.precision = mp.mp.dps
            named_const.base.value = Decimal(str(const_func()))
        elif force:
            if verbose:
                print(f'    Skipping calculation of {const}, too inefficient or no calculation available!')
            if PREFIX_LINK not in const_func.__doc__:
                if verbose:
                    print(f'No backup value exists for {const}! Add one when possible (probably using OEIS)')
            else:
                i = const_func.__doc__.index(PREFIX_LINK)
                url = const_func.__doc__[i + len(PREFIX_LINK) : i + const_func.__doc__[i:].index('\n')]
                url = f'{url}/b{url[url.index(PREFIX_URL) + len(PREFIX_URL) : ]}.txt'
                try:
                    value, precision2 = Universal.read_oeis(urlopen(url))
                    if precision2 < MIN_PRECISION:
                        if verbose:
                            print(f'    OEIS value has too low precision, check back if and when it has at least {MIN_PRECISION} digits.')
                    else:
                        if verbose:
                            print('    OEIS value found, will be used instead')
                        named_const.base.value = value[:16001] # the numeric type is limited to 16383 digits after the decimal point apparently, so for now this sits here
                        named_const.base.precision = min(precision2, 16000)
                except:
                    if verbose:
                        print(f'Exception while fetching {const} from OEIS: {format_exc()}')
        
        if isinstance(name, str):
            named_const.name = const
            named_const.description = const_func.__doc__[:const_func.__doc__.index('.\n')].lstrip()
        if force or named_const.value:
            return named_const
    
    @staticmethod
    def calc_derived(db, ext: Tuple[str, dict] or DerivedConstant, base=None):
        family, args = ext if isinstance(ext, tuple) else ext.family, ext.args
        if family not in DerivedConstants.__dict__:
            return None # who dis
        
        res = DerivedConstants.__dict__[family].__get__(0)(db, **args)
        res = res if isinstance(res, mp.mpf) else res.value
        
        derived = ext if isinstance(ext, DerivedConstant) else DerivedConstant()
        if not derived.base:
            derived.base = base if base else Constant()
        derived.base.value = Decimal(str(res))
        if isinstance(ext, tuple):
            derived.family, derived.args = ext
        return derived
    
    @staticmethod
    def fill_pcf_canonical(const: PcfCanonicalConstant, pcf: PCF, calculation: PCFCalc or None = None):
        const.original_a = [int(coef) for coef in pcf.a.all_coeffs()]
        const.original_b = [int(coef) for coef in pcf.b.all_coeffs()]
        top, bot = pcf.get_canonical_form()
        const.P = [int(coef) for coef in top.all_coeffs()]
        const.Q = [int(coef) for coef in bot.all_coeffs()]
        if calculation:
            getcontext().prec = min(calculation.precision + 10, 16000)
            const.base.value = Decimal(str(calculation.value))
            const.base.precision = calculation.precision
            const.last_matrix = reduce(lambda a, b: a + ',' + str(b), calculation.last_matrix[1:], str(calculation.last_matrix[0]))
            const.depth = calculation.depth
            const.convergence = calculation.convergence
        return const
    
    @staticmethod
    def calc_pcf(ext: Tuple[PCF, List[int] or None, int or None] or PcfCanonicalConstant, base=None, depth_multiplier=None):
        pcf = None
        if isinstance(ext, tuple):
            pcf, prev, depth = ext
            ext = PcfCanonicalConstant()
            ext.base = Constant()
            depth_multiplier = depth_multiplier if depth_multiplier else 1 # when given explicitly, the user will probably want the exact depth they're specifying
        else:
            pcf = PCF(ext.original_a, ext.original_b) if ext.original_a else PCF.from_canonical_form((ext.P, ext.Q))
            prev = [int(x) for x in ext.last_matrix.split(',')]
            depth = ext.depth
            depth_multiplier = depth_multiplier if depth_multiplier else 2 # intended for use from calc or calc_silent, then the multiplier gives this meaning
        
        res = PCFCalc(pcf, prev, depth).run(depth=depth*depth_multiplier, precision=0)
        return Universal.fill_pcf_canonical(ext, pcf, res)
    
    @staticmethod
    def calc_silent(ext, db, const=None):
        if const:
            Universal.set_precision(const.precision)
        
        if isinstance(ext, NamedConstant):
            return Universal.calc_named(ext, const)
        if isinstance(ext, DerivedConstant):
            return Universal.calc_derived(db, ext, const)
        if isinstance(ext, PcfCanonicalConstant):
            return Universal.calc_pcf(ext, const)
        # returning nothing because unrecognized!
    
    @staticmethod
    def calc(ext, db, const=None):
        res = Universal.calc_silent(ext, db, const)
        if not res:
            raise Exception(f'Error while calculating constant with extension of type {type(ext)}: Unrecognized type or calculation failed')
        return res

class DerivedConstants:

    @staticmethod
    def set_precision(prec: int = 4000) -> None:
        '''
        set the precision (in significant digits in base 10). shared with other classes in this file
        '''
        mp.mp.dps = prec
    
    @staticmethod
    def lerch(db, z, s, a):
        '''
        computes the lerch transcendent phi(z,s,a), see https://en.wikipedia.org/wiki/Lerch_zeta_function
        '''
        return mp.lerchphi(z, s, a)
    
    @staticmethod
    def euler_pcf(db, base_roots, diffs, timeout_sec=0):
        '''
        computes the euler PCF defined by base_roots and diffs, with an optional timeout.
        returns a PCFCalc.Result instead of an mp.mpf
        p.s.: zigzagzeta is a special case of this
        '''
        if len(base_roots) != len(diffs):
            raise Exception('base_roots and diffs must have same length')
        n = Symbol('n')
        h1 = reduce(mul, [n - r for r in base_roots])
        h2 = reduce(mul, [n - r - d for r, d in zip(base_roots, diffs)])
        return PCFCalc(PCF(Poly(h1 + compose(h2, n + 1)), Poly(-h1 * h2))).run(precision = mp.mp.dps, timeout_sec = timeout_sec)

    @staticmethod
    def pcf_family_instance(db, family_id, args, timeout_sec=0):
        '''
        evaluates the PCF from the given family, using the given arguments.
        the string representation of the family's polynomials is expected to have its arguments as c0,c1,c2,c3,...
        also not all arguments have to be present in both polynomials in the family.
        returns a PCFCalc.Result instead of an mp.mpf
        '''
        family = db.session.query(PcfFamily).filter(PcfFamily.family_id == family_id).first()
        if not family:
            raise Exception(f'pcf family with id {family_id} not found')
        a = Poly(family.a).eval({f'c{i}': v for i, v in enumerate(args) if f'c{i}' in family.a})
        b = Poly(family.b).eval({f'c{i}': v for i, v in enumerate(args) if f'c{i}' in family.b})
        return PCFCalc(PCF(a, b)).run(precision = mp.mp.dps, timeout_sec = timeout_sec)

    @staticmethod
    def lwt1(db, aux, poly, start):
        '''
        computes the lindemann-weierstrass transcendental using the given auxillary function,
        string representation of a bivariate polynomial poly(x,a), and initial guess for the root of poly(x,aux(x)).
        the result is guaranteed to be transcendental iff the auxillary function is transcendental, and the result is nonzero.
        see https://en.wikipedia.org/wiki/Lindemann%E2%80%93Weierstrass_theorem
        
        raises an exception if aux is not a recognized function
        '''
        if aux not in mp.__dict__:
            raise Exception(f'unrecognized auxillary function {aux}, must be implemented in mpmath')
        poly = Poly(poly)
        aux = mp.__dict__[aux]
        return mp.findroot(lambda x : poly.eval({'x' : x, 'a' : aux(x)}), start)

    @staticmethod
    def gst(db, base: float or str, power: float or str):
        '''
        computes the gelfond-schneider transcendental base ** power. here, base and power can each
        either be numeric values, or can be strings which represent const_ids in the database, which will be automatically queried.
        the result is guaranteed to be transcendental iff base is neither 0 nor 1, and power is irrational.
        see https://en.wikipedia.org/wiki/Gelfond%E2%80%93Schneider_theorem
        
        raises an exception if either base or power are strings which aren't valid const_ids in the database
        '''
        if isinstance(base, str):
            base_const = db.session.query(Constant).filter(Constant.const_id == base).first()
            if not base_const:
                raise Exception(f'constant with id {base} not found')
            base = mp.mpf(str(base_const.value))
        if isinstance(power, str):
            power_const = db.session.query(Constant).filter(Constant.const_id == power).first()
            if not power_const:
                raise Exception(f'constant with id {power} not found')
            power = mp.mpf(str(power_const.value))
        return mp.mpf(base) ** power
    
    @staticmethod
    def mpmath(db, func, **args):
        if func not in mp.__dict__:
            raise Exception(f'Function {func} is not part of the mpmath library')
        return mp.__dict__[func](**args)
        

# TODO use sympy for primes!
# TODO consider using gosper's acceleration of series?
# TODO in general, nprod and nsum seem to be faster when using direct method (and sometimes this is even the only correct way). investigate other methods sometime?
class Constants:
    '''
    arbitrary-precision calculations of named constants.
    
    This class aims to contain most of https://en.wikipedia.org/wiki/List_of_mathematical_constants, 
    excluding rationals, non-reals, and redundant constants (which are connected
    via (2,2)-degree relations to other constants already in here).
    
    PLEASE NOTE: the comments here are formatted in a very specific way,
    meant to serve `create_db.py`! Do not tamper with the comments!
    '''
    
    # If there's a WARNING or a CAUTION, it's a constant that takes a long (or somewhat long)
    # time to calculate to 4000 precision. Find ways to calculate it more efficiently, if possible.
    # If there's a TODO, it's a constant that needs to be added. Do that sometime.

    @staticmethod
    def set_precision(prec: int = 4000) -> None:
        '''
        set the precision (in significant digits in base 10). shared with other classes in this file
        '''
        mp.mp.dps = prec

    @staticmethod
    def pi() -> mpf:
        '''
        pi, fundamental circle constant.
        OEIS link: https://oeis.org/A000796
        '''
        return mp.pi()
    
    @staticmethod
    def sqrt2() -> mpf:
        '''
        square root of 2, also called pythagoras constant.
        OEIS link: https://oeis.org/A002193
        '''
        return mp.sqrt(2)
    
    @staticmethod
    def sqrt3() -> mpf:
        '''
        square root of 3, also called theodorus constant.
        OEIS link: https://oeis.org/A002194
        '''
        return mp.sqrt(3)
    
    @staticmethod
    def phi() -> mpf:
        '''
        golden ratio, positive root of phi^2 - phi - 1.
        OEIS link: https://oeis.org/A001622
        '''
        return mp.phi()
    
    @staticmethod
    def cbrt2() -> mpf:
        '''
        cube root of 2, related to doubling cubes.
        OEIS link: https://oeis.org/A002580
        '''
        return mp.cbrt(2)
    
    @staticmethod
    def cbrt3() -> mpf:
        '''
        cube root of 3.
        OEIS link: https://oeis.org/A002581
        '''
        return mp.cbrt(3)
    
    @staticmethod
    def root12of2() -> mpf:
        '''
        12th root of 2, basis of modern western music theory.
        OEIS link: https://oeis.org/A010774
        '''
        return mp.root(2, 12)
    
    @staticmethod
    def psi() -> mpf:
        '''
        supergolden ratio, real root of psi^3 - psi^2 - 1.
        OEIS link: https://oeis.org/A092526
        '''
        r = 3 * mp.sqrt(93)
        p1 = mp.cbrt((29 + r) / 2)
        p2 = mp.cbrt((29 - r) / 2)
        return (1 + p1 + p2) / 3
    
    @staticmethod
    def mu() -> mpf:
        '''
        hexagonal lattice connective constant, largest root of mu^4 - 4mu^2 + 2.
        OEIS link: https://oeis.org/A179260
        '''
        return mp.sqrt(2 + mp.sqrt(2))
    
    @staticmethod
    def Kprime() -> mpf:
        '''
        kepler bouwkamp constant, also called the polygon inscribing constant.
        OEIS link: https://oeis.org/A085365
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        def iteration(k):
            z = mp.zeta(2 * k)
            pow2 = mp.power(2, 2 * k)
            return (pow2 - 1) / (2 * k) * z * (z - 1 - 1 / pow2)
        return mp.exp(-2 * mp.nsum(iteration, [1, mp.inf], method='d'))
    
    @staticmethod
    def W() -> mpf:
        '''
        wallis constant, real root of w^3 - 2w - 5.
        OEIS link: https://oeis.org/A007493
        '''
        r = mp.sqrt(1929)
        p1 = mp.cbrt((45 + r) / 18)
        p2 = mp.cbrt((45 - r) / 18)
        return p1 + p2
    
    @staticmethod
    def e() -> mpf:
        '''
        euler number, base of the natural logarithm.
        OEIS link: https://oeis.org/A001113
        '''
        return mp.e()
    
    @staticmethod
    def ln2() -> mpf:
        '''
        natural log of 2, has many series representations, and appears often in other constants.
        OEIS link: https://oeis.org/A002162
        '''
        return mp.ln(2)
    
    @staticmethod
    def G_025() -> mpf:
        '''
        gamma(0.25), appears often in other constants.
        OEIS link: https://oeis.org/A068466
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.gamma(0.25)
    
    @staticmethod
    def gamma() -> mpf:
        '''
        euler mascheroni constant, relating the harmonic series and the natural log.
        OEIS link: https://oeis.org/A001620
        '''
        return mp.euler()
    
    @staticmethod
    def E() -> mpf:
        '''
        erdos borwein constant, related to mersenne numbers.
        OEIS link: https://oeis.org/A065442
        '''
        def iteration(n):
            pow2 = mp.power(2, n)
            return mp.power(pow2, -n) * (pow2 + 1) / (pow2 - 1)
        return mp.nsum(iteration, [1, mp.inf], method='d')
    
    @staticmethod
    def Omega() -> mpf:
        '''
        omega constant, real root of omega * e^omega - 1.
        OEIS link: https://oeis.org/A030178
        '''
        return mp.lambertw(1)
    
    @staticmethod
    def Zeta3() -> mpf:
        '''
        apery constant, appears often in physics.
        OEIS link: https://oeis.org/A002117
        '''
        return mp.apery()
    
    @staticmethod
    def L_lim() -> mpf:
        '''
        laplace limit, important to kepler's equation.
        OEIS link: https://oeis.org/A033259
        '''
        def equation(x):
            s = mp.hypot(x, 1)
            return x * mp.exp(s) - s - 1
        return mp.findroot(equation, 0.66)
    
    @staticmethod
    def R_S() -> mpf:
        '''
        ramanujan soldner constant, central to the logarithmic integral.
        OEIS link: https://oeis.org/A070769
        '''
        return mp.findroot(mp.li, 1.5)
    
    @staticmethod
    def G() -> mpf:
        '''
        gauss constant, related to bernoulli's lemniscate.
        OEIS link: https://oeis.org/A014549
        '''
        return 1 / mp.agm(1, mp.sqrt(2))
    
    @staticmethod
    def L_1() -> mpf:
        '''
        first lemniscate constant, related to bernoulli's lemniscate.
        OEIS link: http://oeis.org/A085565
        '''
        return Constants.G() * mp.pi / 2
    
    @staticmethod
    def L_2() -> mpf:
        '''
        second lemniscate constant, related to bernoulli's lemniscate.
        OEIS link: http://oeis.org/A076390
        '''
        return 0.5 / Constants.G()
    
    @staticmethod
    def L() -> mpf:
        '''
        liouville constant, a special case of liouville numbers.
        OEIS link: https://oeis.org/A012245
        '''
        return mp.nsum(lambda n: mp.power(10, -mp.fac(n)), [1, mp.inf], method='d')
    
    @staticmethod
    def C_1() -> mpf:
        '''
        first continued fraction constant.
        OEIS link: https://oeis.org/A052119
        '''
        return mp.besseli(1, 2) / mp.besseli(0, 2)
    
    @staticmethod
    def R() -> mpf:
        '''
        ramanujan constant, infamous almost-integer.
        OEIS link: https://oeis.org/A060295
        '''
        return mp.exp(mp.pi * mp.sqrt(163))
    
    @staticmethod
    def A() -> mpf:
        '''
        glaisher kinkelin constant, related to gamma functions and zeta functions.
        OEIS link: https://oeis.org/A074962
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.glaisher()
    
    @staticmethod
    def C() -> mpf:
        '''
        catalan constant, important to combinatorics, topology, and more.
        OEIS link: https://oeis.org/A006752
        '''
        return mp.catalan()
    
    @staticmethod
    def D() -> mpf:
        '''
        dottie number, real root of cos(d) - d (in radians).
        OEIS link: https://oeis.org/A003957
        '''
        return mp.findroot(lambda x: mp.cos(x) - x, 0.74)
    
    @staticmethod
    def M() -> mpf:
        '''
        meissel mertens constant, one of many constants relating prime numbers.
        OEIS link: https://oeis.org/A077761
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.mertens()
    
    @staticmethod
    def P() -> mpf:
        '''
        universal parabolic constant, a fundamental ratio of parabolas.
        OEIS link: https://oeis.org/A103710
        '''
        sqrt2 = mp.sqrt(2)
        return mp.ln(1 + sqrt2) + sqrt2
    
    @staticmethod
    def C_Cahen() -> mpf:
        '''
        cahen constant, related to the sylvester sequence.
        OEIS link: https://oeis.org/A118227
        '''
        sylvester_dict = dict() # caching the sylvester sequence makes this way faster
        def sylvester(k):
            if k in sylvester_dict:
                return sylvester_dict[k]
            res = mpf(2) if k == 0 else 1 + mp.nprod(sylvester, [0, k - 1])
            sylvester_dict[k] = res
            return res
        return mp.nsum(lambda k: (-1)**k / (sylvester(k) - 1), [0, mp.inf], method='d')
    
    @staticmethod
    def epi() -> mpf:
        '''
        gelfond constant, a result of the gelfond-schneider theorem.
        OEIS link: https://oeis.org/A039661
        '''
        return mp.exp(mp.pi)
    
    @staticmethod
    def G_S() -> mpf:
        '''
        gelfond schneider constant, also called the hilbert number.
        OEIS link: https://oeis.org/A007507
        '''
        return mp.power(2, mp.sqrt(2))
    
    @staticmethod
    def g() -> mpf:
        '''
        golden angle, related to the golden ratio.
        OEIS link: https://oeis.org/A131988
        '''
        return 2 * mp.pi / (1 + mp.phi)
    
    @staticmethod
    def S() -> mpf:
        '''
        sierpinski constant, related to gauss constant and euler-mascheroni constant.
        OEIS link: https://oeis.org/A062089
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.pi * (2 * mp.ln(2) + 3 * mp.ln(mp.pi) + 2 * mp.euler - 4 * mp.ln(mp.gamma(0.25)))
    
    @staticmethod
    def L_R() -> mpf:
        '''
        landau ramanujan constant, central to a theorem by edmund landau.
        OEIS link: https://oeis.org/A064533
        WARNING: This is not a calculation!
        ''' # TODO implement prime sieve, then iterate over primes congruent 1 mod 4...
        return mp.mpf('0.76422365358922066299')
    
    @staticmethod
    def G_L() -> mpf:
        '''
        gieseking constant, also called lobachevsky constant.
        OEIS link: https://oeis.org/A143298
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.clsin(2 , mp.pi / 3)
    
    @staticmethod
    def beta() -> mpf:
        '''
        bernstein constant, describing errors of best uniform approximations.
        OEIS link: https://oeis.org/A073001
        WARNING: This is not a calculation!
        ''' # TODO implement Varga&Carpenter algorithm?
        return mp.mpf('0.28016949902386913303')
    
    @staticmethod
    def T() -> mpf:
        '''
        tribonacci constant, real root of t^3 - t^2 - t - 1.
        OEIS link: https://oeis.org/A058265
        '''
        r = 3 * mp.sqrt(33)
        p1 = mp.cbrt(19 + r)
        p2 = mp.cbrt(19 - r)
        return (1 + p1 + p2) / 3
    
    @staticmethod
    def B_2() -> mpf:
        '''
        brun constant, follows from brun's theorem.
        OEIS link: https://oeis.org/A065421
        WARNING: This is not a calculation!
        ''' # TODO need primes again... might not include this after all, since to get 13 significant digits you need all twin primes up to 10^16!
        return mp.mpf('1.902160583104')
    
    @staticmethod
    def Pi_2() -> mpf:
        '''
        twin primes constant, central to the twin primes conjecture.
        OEIS link: https://oeis.org/A005597
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.twinprime()
    
    @staticmethod
    def rho() -> mpf:
        '''
        plastic number, real root of rho^3 - rho - 1.
        OEIS link: https://oeis.org/A060006
        '''
        r = mp.sqrt(69) / 18
        p1 = mp.cbrt(0.5 + r)
        p2 = mp.cbrt(0.5 - r)
        return p1 + p2
    
    @staticmethod
    def z_975() -> mpf:
        '''
        z score for 97.5 percentile point, commonly used alongside normal distributions.
        OEIS link: https://oeis.org/A220510
        '''
        return mp.sqrt(2) * mp.erfinv(0.95)
    
    @staticmethod
    def tau() -> mpf:
        '''
        prouhet thue morse constant, appears in probability.
        OEIS link: https://oeis.org/A014571
        '''
        return 0.25 * (2 - mp.nprod(lambda n: 1 - mp.power(2, -mp.power(2, n)), [0, mp.inf]))
    
    @staticmethod
    def lambda_GD() -> mpf:
        '''
        golomb dickman constant, appears in random permutation theory and number theory.
        OEIS link: https://oeis.org/A084945
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.quad(lambda x: mp.exp(mp.li(x)), [0, 1])
    
    @staticmethod
    def c() -> mpf:
        '''
        asymptotic lebesgue constant.
        OEIS link: https://oeis.org/A243277
        WARNING: This is not a calculation!
        '''
        ## TODO this code seems right but gives wrong result???
        #s = -mp.digamma(0.5)
        #return 4 / mp.pi ** 2 * (mp.nsum(lambda k: 2 * mp.ln(k) / (4 * mp.power(k, 2) - 1), [1, mp.inf], method='d') + s)
        return mp.mpf('0.98943127383114695174')
    
    @staticmethod
    def C_FT() -> mpf:
        '''
        feller tornier constant, describing certain prime factorizations.
        OEIS link: https://oeis.org/A065493
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.66131704946962233528')
    
    @staticmethod
    def C_10() -> mpf:
        '''
        base10 champernowne constant.
        OEIS link: https://oeis.org/A033307
        '''
        res = '0.'
        i = 1
        while len(res) < mp.mp.dps:
            res += str(i)
            i += 1
        return mpf(res)
    
    @staticmethod
    def sigma_10() -> mpf:
        '''
        salem constant, smallest known salem number.
        OEIS link: https://oeis.org/A073011
        '''
        return mp.findroot(lambda x: mp.polyval([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], x), 1.2)
    
    @staticmethod
    def K_0() -> mpf:
        '''
        khinchin constant, a surprising fundamental constant in continued fractions.
        OEIS link: https://oeis.org/A002210
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.khinchin()
    
    @staticmethod
    def beta_Levy() -> mpf:
        '''
        first levy constant, related to asymptotic behavior in continued fractions.
        OEIS link: https://oeis.org/A100199
        '''
        return mp.pi ** 2 / (12 * mp.ln(2))
    
    @staticmethod
    def eLevy() -> mpf:
        '''
        second levy constant, related to asymptotic behavior in continued fractions.
        OEIS link: https://oeis.org/A086702
        '''
        return mp.exp(Constants.beta_Levy())
    
    @staticmethod
    def C_CE() -> mpf:
        '''
        copeland erdos constant.
        OEIS link: https://oeis.org/A033308
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.23571113171923293137')
    
    @staticmethod
    def A_Pi() -> mpf:
        '''
        mills constant, "smallest" real that generates prime numbers via exponents.
        OEIS link: https://oeis.org/A051021
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated
        return mp.mpf('1.30637788386308069046')
    
    @staticmethod
    def delta_G() -> mpf:
        '''
        gompertz constant, appears in some special integrals.
        OEIS link: https://oeis.org/A073003
        '''
        return -mp.e * mp.ei(-1)
    
    @staticmethod
    def V_dp() -> mpf:
        '''
        van der pauw constant, involved in the van der pauw method.
        OEIS link: https://oeis.org/A163973
        '''
        return mp.pi / mp.ln(2)
    
    @staticmethod
    def theta_m() -> mpf:
        '''
        magic angle, important to magnetic resonance imaging.
        OEIS link: https://oeis.org/A195696
        '''
        return mp.atan(mp.sqrt(2))
    
    @staticmethod
    def C_Artin() -> mpf:
        '''
        artin constant, related to artin's conjecture on primitive roots.
        OEIS link: https://oeis.org/A005596
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.37395581361920228805')
    
    @staticmethod
    def C_P() -> mpf:
        '''
        porter constant, related to the efficiency of euclid algorithm.
        OEIS link: https://oeis.org/A086237
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        ln2 = mp.ln(2)
        pi2 = mp.pi ** 2
        return 6 * ln2 / pi2 * (3 * ln2 + 4 * mp.euler - 24 / pi2 * mp.zeta(2, derivative=1) - 2) - 0.5
    
    @staticmethod
    def L_Lochs() -> mpf:
        '''
        lochs constant, involved in lochs' theorem regarding continued fractions.
        OEIS link: https://oeis.org/A086819
        '''
        return 6 * mp.ln(2) * mp.ln(10) / mp.pi ** 2
    
    @staticmethod
    def D_V() -> mpf:
        '''
        devicci tesseract constant, describing the largest cube that can pass through a 4d hypercube.
        OEIS link: https://oeis.org/A243309
        '''
        return mp.findroot(lambda x: mp.polyval([4, 0, -28, 0, -7, 0, 16, 0, 16], x), 1)
    
    @staticmethod
    def C_N() -> mpf:
        '''
        niven constant, largest exponent in prime factorizations "on average".
        OEIS link: https://oeis.org/A033150
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return 1 + mp.nsum(lambda n: 1 - 1 / mp.zeta(n), [2, mp.inf], method='d')
    
    @staticmethod
    def S_Pi() -> mpf:
        '''
        stephens constant, density of some subsets of primes.
        OEIS link: https://oeis.org/A065478
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.57595996889294543964')
    
    @staticmethod
    def P_Dragon() -> mpf:
        '''
        paperfolding constant, related to the dragon curve.
        OEIS link: https://oeis.org/A143347
        '''
        def iteration(n):
            two_n = mp.power(2, n)
            return mp.power(2, -two_n) / (1 - mp.power(2, -4 * two_n))
        return mp.nsum(iteration, [0, mp.inf], method='d')
    
    @staticmethod
    def psi_Fib() -> mpf:
        '''
        reciprocal fibonacci constant, sum of reciprocals of fibonacci numbers.
        OEIS link: https://oeis.org/A079586
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.nsum(lambda n: 1 / mp.fib(n), [1, mp.inf], method='d')
    
    @staticmethod
    def delta() -> mpf:
        '''
        first feigenbaum constant, important to bifurcation theory.
        OEIS link: https://oeis.org/A006890
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated? maybe https://rosettacode.org/wiki/Feigenbaum_constant_calculation#Python
        return mp.mpf('4.66920160910299067185')
    
    @staticmethod
    def Delta_3() -> mpf:
        '''
        robbins constant, mean length of random line segments in a unit cube.
        OEIS link: https://oeis.org/A073012
        '''
        sqrt2 = mp.sqrt(2)
        sqrt3 = mp.sqrt(3)
        p1 = (4 + 17 * sqrt2 - 6 * sqrt3 - 7 * mp.pi) / 105
        p2 = mp.ln(1 + sqrt2) / 5
        p3 = mp.ln(2 + sqrt3) * 2 / 5
        return p1 + p2 + p3
    
    @staticmethod
    def W_S() -> mpf:
        '''
        weierstrass constant.
        OEIS link: https://oeis.org/A094692
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.power(2, 1.25) * mp.sqrt(mp.pi) * mp.exp(mp.pi / 8) / mp.gamma(0.25) ** 2
    
    @staticmethod
    def F() -> mpf:
        '''
        fransen robinson constant, related to the reciprocal gamma function.
        OEIS link: https://oeis.org/A058655
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.quad(mp.rgamma, [0, mp.inf])
    
    @staticmethod
    def alpha() -> mpf:
        '''
        second feigenbaum contsant, important to bifurcation theory.
        OEIS link: https://oeis.org/A006891
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated? maybe https://github.com/brorson/FeigenbaumConstants
        return mp.mpf('2.50290787509589282228')
    
    @staticmethod
    def C_2() -> mpf:
        '''
        second du bois reymond constant.
        OEIS link: https://oeis.org/A062546
        '''
        return (mp.exp(2) - 7) / 2
    
    @staticmethod
    def delta_ETF() -> mpf:
        '''
        erdos tenenbaum ford constant, appears in number theory.
        OEIS link: https://oeis.org/A074738
        '''
        ln2 = mp.ln(2)
        return 1 - (1 + mp.ln(ln2)) / ln2
    
    @staticmethod
    def lambda_C() -> mpf:
        '''
        conway constant, related to the look-and-say sequence.
        OEIS link: https://oeis.org/A014715
        '''
        return mp.findroot(lambda x: mp.polyval([1, 0, -1, -2, -1, 2, 2, 1, -1, -1, -1, -1, -1, 2, 5, 3, -2, -10, -3, -2, 6, 6, 1, 9, -3, -7, -8, -8, 10, 6, 8, -5, -12, 7, -7, 7, 1, -3, 10, 1, -6, -2, -10, -3, 2, 9, -3, 14, -8, 0, -7, 9, 3, -4, -10, -7, 12, 7, 2, -12, -4, -2, 5, 0, 1, -7, 7, -4, 12, -6, 3, -6], x) , 1.3)
    
    @staticmethod
    def sigma() -> mpf:
        '''
        hafner sarnak mccurley constant, related to coprime determinants of integer matrices.
        OEIS link: https://oeis.org/A085849
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.35323637185499598454')
    
    @staticmethod
    def B_H() -> mpf:
        '''
        backhouse constant, constructed using power series with prime coefficients.
        OEIS link: https://oeis.org/A072508
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated
        return mp.mpf('1.45607494858268967139')
    
    @staticmethod
    def V() -> mpf:
        '''
        viswanath constant, related to random fibonacci sequences.
        OEIS link: https://oeis.org/A078416
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated
        return mp.mpf('1.1319882487943')
    
    @staticmethod
    def q() -> mpf:
        '''
        komornik loreti constant, related to non-integer representations.
        OEIS link: https://oeis.org/A055060
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.findroot(lambda q: mp.nprod(lambda n: 1 - mp.power(q, -mp.power(2, n)), [0, mp.inf]) + (q - 2) / (q - 1), 2)
    
    @staticmethod
    def C_HBM() -> mpf:
        '''
        heath brown moroz constant, related to the cubic surface w^3 = xyz.
        OEIS link: https://oeis.org/A118228
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.00131764115485317810')
    
    @staticmethod
    def S_MRB() -> mpf:
        '''
        mrb constant, named after marvin ray burns.
        OEIS link: https://oeis.org/A037077
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.nsum(lambda n: mp.power(-1, n) * (mp.root(n, n) - 1), [1, mp.inf])
    
    @staticmethod
    def rho_Pi() -> mpf:
        '''
        prime constant, constructed from indicators of prime numbers.
        OEIS link: https://oeis.org/A051006
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.41468250985111166024')
    
    @staticmethod
    def sigma_S() -> mpf:
        '''
        somos quadratic recurrence constant, related to the lerch transcendent.
        OEIS link: https://oeis.org/A112302
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.nprod(lambda n: mp.power(n, mp.power(2, -n)), [1, mp.inf], method='d')
    
    @staticmethod
    def alpha_F() -> mpf:
        '''
        foias constant, only number for which a certain recurrence diverges.
        OEIS link: https://oeis.org/A085848
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated (maybe findroot???)
        return mp.mpf('1.18745235112650105459')
    
    @staticmethod
    def L_D() -> mpf:
        '''
        unit disk logarithmic capacity.
        OEIS link: https://oeis.org/A249205
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.power(mp.gamma(0.25), 2) / (4 * mp.power(mp.pi, 1.5))
    
    @staticmethod
    def T_Pi() -> mpf:
        '''
        taniguchi constant, a kind of euler product.
        OEIS link: https://oeis.org/A175639
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.67823449191739197803')

'''
List of rejected constants:
    No rationals allowed:
        1, 2, 0.5, 0, -1
    No complex numbers allowed:
        i
    No known method to compute:
        Bloch constant, Landau constant, Third Landau constant, de Bruijn-Newman contant,
        Binary Alphabet Chvátal-Sankoff constant, Chaitin constant(s) (in general) (sorta),
        Embree-Trefethen constant
    Related to constants already added:
        sqrt(5): (1,1)-degree relation with phi
        Silver ratio: (1,1)-degree relation with sqrt(2)
        Lemniscate constant: (1,1)-degree relation with First Lemniscate constant
        Second Hermite constant: (2,1)-degree relation with sqrt(3)
        Second Favard constant: (2,2)-degree relation with pi
        First Nielsen-Ramanujan constant: (2,2)-degree relation with pi
        Lieb square ice constant: (2,1)-degree relation with sqrt(3)
'''
