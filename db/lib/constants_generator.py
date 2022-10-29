import numpy as np
import mpmath as mp
from mpmath import mpf
from decimal import Decimal, getcontext
from typing import List

# TODO use sympy for primes!
# TODO consider using gosper's acceleration of series?
# TODO in general, nprod and nsum seem to be faster when using direct method (and sometimes this is even the only correct way). investigate other methods sometime?
class Constants:
    '''
    arbitrary-precision calculations of constants.
    
    This class aims to contain most of https://en.wikipedia.org/wiki/List_of_mathematical_constants, 
    excluding rationals, non-reals, and redundant constants (which are connected
    via (2,2)-degree relations to other constants already in here).
    '''
    
    # If there's a WARNING or a CAUTION, it's a constant that takes a long (or somewhat long)
    # time to calculate to 4000 precision. Find ways to calculate it more efficiently, if possible.
    # If there's a TODO, it's a constant that needs to be added. Do that sometime.

    @staticmethod
    def set_precision(prec: int = 4000) -> None:
        '''
        set the precision (in significant digits in base 10).
        '''
        getcontext().prec = prec # might be redundant...
        mp.mp.dps = prec

    @staticmethod
    def pi() -> mpf:
        '''
        pi, fundamental circle constant.
        '''
        return mp.pi()
    
    @staticmethod
    def sqrt2() -> mpf:
        '''
        square root of 2, also called pythagoras constant.
        '''
        return mp.sqrt(2)
    
    @staticmethod
    def sqrt3() -> mpf:
        '''
        square root of 3, also called theodorus constant.
        '''
        return mp.sqrt(3)
    
    @staticmethod
    def phi() -> mpf:
        '''
        golden ratio, positive root of phi^2 - phi - 1.
        '''
        return mp.phi()
    
    @staticmethod
    def cbrt2() -> mpf:
        '''
        cube root of 2, related to doubling cubes.
        '''
        return mp.cbrt(2)
    
    @staticmethod
    def cbrt3() -> mpf:
        '''
        cube root of 3.
        '''
        return mp.cbrt(3)
    
    @staticmethod
    def root12of2() -> mpf:
        '''
        12th root of 2, basis of modern western music theory.
        '''
        return mp.root(2, 12)
    
    @staticmethod
    def psi() -> mpf:
        '''
        supergolden ratio, real root of psi^3 - psi^2 - 1.
        '''
        r = 3 * mp.sqrt(93)
        p1 = mp.cbrt((29 + r) / 2)
        p2 = mp.cbrt((29 - r) / 2)
        return (1 + p1 + p2) / 3
    
    @staticmethod
    def mu() -> mpf:
        '''
        hexagonal lattice connective constant, largest root of mu^4 - 4mu^2 + 2.
        '''
        return mp.sqrt(2 + mp.sqrt(2))
    
    @staticmethod
    def Kprime() -> mpf:
        '''
        kepler bouwkamp constant, also called the polygon inscribing constant.
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
        '''
        r = mp.sqrt(1929)
        p1 = mp.cbrt((45 + r) / 18)
        p2 = mp.cbrt((45 - r) / 18)
        return p1 + p2
    
    @staticmethod
    def e() -> mpf:
        '''
        euler number, base of the natural logarithm.
        '''
        return mp.e()
    
    @staticmethod
    def ln2() -> mpf:
        '''
        natural log of 2, has many series representations, and appears often in other constants.
        '''
        return mp.ln(2)
    
    @staticmethod
    def G025() -> mpf:
        '''
        gamma(0.25), appears often in other constants.
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.gamma(0.25)
    
    @staticmethod
    def gamma() -> mpf:
        '''
        euler mascheroni constant, relating the harmonic series and the natural log.
        '''
        return mp.euler()
    
    @staticmethod
    def E() -> mpf:
        '''
        erdos borwein constant, related to mersenne numbers.
        '''
        def iteration(n):
            pow2 = mp.power(2, n)
            return mp.power(pow2, -n) * (pow2 + 1) / (pow2 - 1)
        return mp.nsum(iteration, [1, mp.inf], method='d')
    
    @staticmethod
    def Omega() -> mpf:
        '''
        omega constant, real root of omega * e^omega - 1.
        '''
        return mp.lambertw(1)
    
    @staticmethod
    def Zeta3() -> mpf:
        '''
        apery constant, appears often in physics.
        '''
        return mp.apery()
    
    @staticmethod
    def L_lim() -> mpf:
        '''
        laplace limit, important to kepler's equation.
        '''
        def equation(x):
            s = mp.hypot(x, 1)
            return x * mp.exp(s) - s - 1
        return mp.findroot(equation, 0.66)
    
    @staticmethod
    def R_S() -> mpf:
        '''
        ramanujan soldner constant, central to the logarithmic integral.
        '''
        return mp.findroot(mp.li, 1.5)
    
    @staticmethod
    def G() -> mpf:
        '''
        gauss constant, related to bernoulli's lemniscate.
        '''
        return 1 / mp.agm(1, mp.sqrt(2))
    
    @staticmethod
    def L1() -> mpf:
        '''
        first lemniscate constant, related to bernoulli's lemniscate.
        '''
        return Constants.G() * mp.pi / 2
    
    @staticmethod
    def L2() -> mpf:
        '''
        second lemniscate constant, related to bernoulli's lemniscate.
        '''
        return 0.5 / Constants.G()
    
    @staticmethod
    def L() -> mpf:
        '''
        liouville constant, a special case of liouville numbers.
        '''
        return mp.nsum(lambda n: mp.power(10, -mp.fac(n)), [1, mp.inf], method='d')
    
    @staticmethod
    def C_1() -> mpf:
        '''
        first continued fraction constant.
        '''
        return mp.besseli(1, 2) / mp.besseli(0, 2)
    
    @staticmethod
    def R() -> mpf:
        '''
        ramanujan constant, infamous almost-integer.
        '''
        return mp.exp(mp.pi * mp.sqrt(163))
    
    @staticmethod
    def A() -> mpf:
        '''
        glaisher kinkelin constant, related to gamma functions and zeta functions.
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.glaisher()
    
    @staticmethod
    def C() -> mpf:
        '''
        catalan constant, important to combinatorics, topology, and more.
        '''
        return mp.catalan()
    
    @staticmethod
    def D() -> mpf:
        '''
        dottie number, real root of cos(d) - d (in radians).
        '''
        return mp.findroot(lambda x: mp.cos(x) - x, 0.74)
    
    @staticmethod
    def M() -> mpf:
        '''
        meissel mertens constant, one of many constants relating prime numbers.
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.mertens()
    
    @staticmethod
    def P() -> mpf:
        '''
        universal parabolic constant, a fundamental ratio of parabolas.
        '''
        sqrt2 = mp.sqrt(2)
        return mp.ln(1 + sqrt2) + sqrt2
    
    @staticmethod
    def C_Cahen() -> mpf:
        '''
        cahen constant, related to the sylvester sequence.
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
        '''
        return mp.exp(mp.pi)
    
    @staticmethod
    def G_S() -> mpf:
        '''
        gelfond schneider constant, also called the hilbert number.
        '''
        return mp.power(2, mp.sqrt(2))
    
    @staticmethod
    def g() -> mpf:
        '''
        golden angle, related to the golden ratio.
        '''
        return 2 * mp.pi / (1 + mp.phi)
    
    @staticmethod
    def S() -> mpf:
        '''
        sierpinski constant, related to gauss constant and euler-mascheroni constant.
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.pi * (2 * mp.ln(2) + 3 * mp.ln(mp.pi) + 2 * mp.euler - 4 * mp.ln(mp.gamma(0.25)))
    
    @staticmethod
    def L_R() -> mpf:
        '''
        landau ramanujan constant, central to a theorem by edmund landau.
        WARNING: This is not a calculation!
        ''' # TODO implement prime sieve, then iterate over primes congruent 1 mod 4...
        return mp.mpf('0.76422365358922066299')
    
    @staticmethod
    def G_L() -> mpf:
        '''
        gieseking constant, also called lobachevsky constant.
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.clsin(2 , mp.pi / 3)
    
    @staticmethod
    def beta() -> mpf:
        '''
        bernstein constant, describing errors of best uniform approximations.
        WARNING: This is not a calculation!
        ''' # TODO implement Varga&Carpenter algorithm?
        return mp.mpf('0.28016949902386913303')
    
    @staticmethod
    def T() -> mpf:
        '''
        tribonacci constant, real root of t^3 - t^2 - t - 1.
        '''
        r = 3 * mp.sqrt(33)
        p1 = mp.cbrt(19 + r)
        p2 = mp.cbrt(19 - r)
        return (1 + p1 + p2) / 3
    
    @staticmethod
    def B_2() -> mpf:
        '''
        brun constant, follows from brun's theorem.
        WARNING: This is not a calculation!
        ''' # TODO need primes again... might not include this after all, since to get 13 significant digits you need all twin primes up to 10^16!
        return mp.mpf('1.902160583104')
    
    @staticmethod
    def Pi2() -> mpf:
        '''
        twin primes constant, central to the twin primes conjecture.
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.twinprime()
    
    @staticmethod
    def rho() -> mpf:
        '''
        plastic number, real root of rho^3 - rho - 1.
        '''
        r = mp.sqrt(69) / 18
        p1 = mp.cbrt(0.5 + r)
        p2 = mp.cbrt(0.5 - r)
        return p1 + p2
    
    @staticmethod
    def z975() -> mpf:
        '''
        z score for 97.5 percentile point, commonly used alongside normal distributions.
        '''
        return mp.sqrt(2) * mp.erfinv(0.95)
    
    @staticmethod
    def tau() -> mpf:
        '''
        prouhet thue morse constant, appears in probability.
        '''
        return 0.25 * (2 - mp.nprod(lambda n: 1 - mp.power(2, -mp.power(2, n)), [0, mp.inf]))
    
    @staticmethod
    def lambda_GD() -> mpf:
        '''
        golomb dickman constant, appears in random permutation theory and number theory.
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.quad(lambda x: mp.exp(mp.li(x)), [0, 1])
    
    @staticmethod
    def c() -> mpf:
        '''
        asymptotic lebesgue constant.
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
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.66131704946962233528')
    
    @staticmethod
    def C_10() -> mpf:
        '''
        base10 champernowne constant.
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
        '''
        return mp.findroot(lambda x: mp.polyval([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], x), 1.2)
    
    @staticmethod
    def K0() -> mpf:
        '''
        khinchin constant, a surprising fundamental constant in continued fractions.
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.khinchin()
    
    @staticmethod
    def beta_Levy() -> mpf:
        '''
        first levy constant, related to asymptotic behavior in continued fractions.
        '''
        return mp.pi ** 2 / (12 * mp.ln(2))
    
    @staticmethod
    def eLevy() -> mpf:
        '''
        second levy constant, related to asymptotic behavior in continued fractions.
        '''
        return mp.exp(Constants.beta_Levy())
    
    @staticmethod
    def C_CE() -> mpf:
        '''
        copeland erdos constant.
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.23571113171923293137')
    
    @staticmethod
    def A_Pi() -> mpf:
        '''
        mills constant, "smallest" real that generates prime numbers via exponents.
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated
        return mp.mpf('1.30637788386308069046')
    
    @staticmethod
    def delta_G() -> mpf:
        '''
        gompertz constant, appears in some special integrals.
        '''
        return -mp.e * mp.ei(-1)
    
    @staticmethod
    def V_dp() -> mpf:
        '''
        van der pauw constant, involved in the van der pauw method.
        '''
        return mp.pi / mp.ln(2)
    
    @staticmethod
    def theta_m() -> mpf:
        '''
        magic angle, important to magnetic resonance imaging.
        '''
        return mp.atan(mp.sqrt(2))
    
    @staticmethod
    def C_Artin() -> mpf:
        '''
        artin constant, related to artin's conjecture on primitive roots.
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.37395581361920228805')
    
    @staticmethod
    def C_P() -> mpf:
        '''
        porter constant, related to the efficiency of euclid algorithm.
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        ln2 = mp.ln(2)
        pi2 = mp.pi ** 2
        return 6 * ln2 / pi2 * (3 * ln2 + 4 * mp.euler - 24 / pi2 * mp.zeta(2, derivative=1) - 2) - 0.5
    
    @staticmethod
    def L_Lochs() -> mpf:
        '''
        lochs constant, involved in lochs' theorem regarding continued fractions.
        '''
        return 6 * mp.ln(2) * mp.ln(10) / mp.pi ** 2
    
    @staticmethod
    def D_V() -> mpf:
        '''
        devicci tesseract constant, describing the largest cube that can pass through a 4d hypercube.
        '''
        return mp.findroot(lambda x: mp.polyval([4, 0, -28, 0, -7, 0, 16, 0, 16], x), 1)
    
    @staticmethod
    def C_N() -> mpf:
        '''
        niven constant, largest exponent in prime factorizations "on average".
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return 1 + mp.nsum(lambda n: 1 - 1 / mp.zeta(n), [2, mp.inf], method='d')
    
    @staticmethod
    def S_Pi() -> mpf:
        '''
        stephens constant, density of some subsets of primes.
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.57595996889294543964')
    
    @staticmethod
    def P_Dragon() -> mpf:
        '''
        paperfolding constant, related to the dragon curve.
        '''
        def iteration(n):
            two_n = mp.power(2, n)
            return mp.power(2, -two_n) / (1 - mp.power(2, -4 * two_n))
        return mp.nsum(iteration, [0, mp.inf], method='d')
    
    @staticmethod
    def psi_Fib() -> mpf:
        '''
        reciprocal fibonacci constant, sum of reciprocals of fibonacci numbers.
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.nsum(lambda n: 1 / mp.fib(n), [1, mp.inf], method='d')
    
    @staticmethod
    def delta() -> mpf:
        '''
        first feigenbaum constant, important to bifurcation theory.
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated? maybe https://rosettacode.org/wiki/Feigenbaum_constant_calculation#Python
        return mp.mpf('4.66920160910299067185')
    
    @staticmethod
    def Delta3() -> mpf:
        '''
        robbins constant, mean length of random line segments in a unit cube.
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
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.power(2, 1.25) * mp.sqrt(mp.pi) * mp.exp(mp.pi / 8) / mp.gamma(0.25) ** 2
    
    @staticmethod
    def F() -> mpf:
        '''
        fransen robinson constant, related to the reciprocal gamma function.
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.quad(mp.rgamma, [0, mp.inf])
    
    @staticmethod
    def alpha() -> mpf:
        '''
        second feigenbaum contsant, important to bifurcation theory.
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated? maybe https://github.com/brorson/FeigenbaumConstants
        return mp.mpf('2.50290787509589282228')
    
    @staticmethod
    def C_2() -> mpf:
        '''
        second du bois reymond constant.
        '''
        return (mp.exp(2) - 7) / 2
    
    @staticmethod
    def delta_ETF() -> mpf:
        '''
        erdos tenenbaum ford constant, appears in number theory.
        '''
        ln2 = mp.ln(2)
        return 1 - (1 + mp.ln(ln2)) / ln2
    
    @staticmethod
    def lambda_C() -> mpf:
        '''
        conway constant, related to the look-and-say sequence.
        '''
        return mp.findroot(lambda x: mp.polyval([1, 0, -1, -2, -1, 2, 2, 1, -1, -1, -1, -1, -1, 2, 5, 3, -2, -10, -3, -2, 6, 6, 1, 9, -3, -7, -8, -8, 10, 6, 8, -5, -12, 7, -7, 7, 1, -3, 10, 1, -6, -2, -10, -3, 2, 9, -3, 14, -8, 0, -7, 9, 3, -4, -10, -7, 12, 7, 2, -12, -4, -2, 5, 0, 1, -7, 7, -4, 12, -6, 3, -6], x) , 1.3)
    
    @staticmethod
    def sigma() -> mpf:
        '''
        hafner sarnak mccurley constant, related to coprime determinants of integer matrices.
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.35323637185499598454')
    
    @staticmethod
    def B_H() -> mpf:
        '''
        backhouse constant, constructed using power series with prime coefficients.
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated
        return mp.mpf('1.45607494858268967139')
    
    @staticmethod
    def V() -> mpf:
        '''
        viswanath constant, related to random fibonacci sequences.
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated
        return mp.mpf('1.1319882487943')
    
    @staticmethod
    def q() -> mpf:
        '''
        komornik loreti constant, related to non-integer representations.
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.findroot(lambda q: mp.nprod(lambda n: 1 - mp.power(q, -mp.power(2, n)), [0, mp.inf]) + (q - 2) / (q - 1), 2)
    
    @staticmethod
    def C_HBM() -> mpf:
        '''
        heath brown moroz constant, related to the cubic surface w^3 = xyz.
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.00131764115485317810')
    
    @staticmethod
    def S_MRB() -> mpf:
        '''
        mrb constant, named after marvin ray burns.
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.nsum(lambda n: mp.power(-1, n) * (mp.root(n, n) - 1), [1, mp.inf])
    
    @staticmethod
    def rho_Pi() -> mpf:
        '''
        prime constant, constructed from indicators of prime numbers.
        WARNING: This is not a calculation!
        ''' # TODO need primes again...
        return mp.mpf('0.41468250985111166024')
    
    @staticmethod
    def sigma_S() -> mpf:
        '''
        somos quadratic recurrence constant, related to the lerch transcendent.
        WARNING: Very inefficient to calculate! Handle with care!
        '''
        return mp.nprod(lambda n: mp.power(n, mp.power(2, -n)), [1, mp.inf], method='d')
    
    @staticmethod
    def alpha_F() -> mpf:
        '''
        foias constant, only number for which a certain recurrence diverges.
        WARNING: This is not a calculation!
        ''' # TODO how the hell is this calculated (maybe findroot???)
        return mp.mpf('1.18745235112650105459')
    
    @staticmethod
    def L_D() -> mpf:
        '''
        unit disk logarithmic capacity.
        CAUTION: Inefficient to calculate, but not as much as the others.
        '''
        return mp.power(mp.gamma(0.25), 2) / (4 * mp.power(mp.pi, 1.5))
    
    @staticmethod
    def T_Pi() -> mpf:
        '''
        taniguchi constant, a kind of euler product.
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
