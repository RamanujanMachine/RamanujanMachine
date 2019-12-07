import numpy as np
from math import gcd, floor
from mpmath import mpf as dec
import mpmath
from sympy import symbols, pprint


class MobiusTransform(object):
    def __init__(self, arr=np.eye(2, dtype=object)):
        """
            This class represents a Mobius Transform stored as:
            ( a b )
            ( c d )
            matrix.
            applying this transform onto x will result in:
            (ax + b) / (cx + d)
        """
        super().__init__()
        self.data = arr

    def __str__(self) -> str:
        return str(self.data)

    def __mul__(self, other):
        """
        Function composition operator.
        but in fact - regular matrix multiplication
        :param other: MobiusTransform object to multiply with (multiply on right)
        :return: result MobiusTransform
        """
        tmp = MobiusTransform(np.matmul(self.data, other.data))
        tmp.normalize()
        return tmp

    def __imul__(self, other):
        """
        self = [self] o [other] (composition
        :param other: mobius transform
        :return: self
        """
        self.data = np.matmul(self.data, other.data)
        self.normalize()
        return self

    def __call__(self, x=None):
        """
        apply transformation on X.
        :param x (for our project it will always be 1):
        :return: result of mobius transformation
        """
        a, b, c, d = self.__values()
        if x is not None:
            numerator = dec(a*x + b)
            denominator = dec(c*x + d)
        else:
            numerator = dec(b)
            denominator = dec(d)
        return numerator / denominator

    def normalize(self):
        """
        to be called after every operation. this should prevent coefficients from exploding.
        """
        a, b, c, d = self.__values()
        divider = gcd(gcd(a, b), gcd(c, d))
        if divider != 1:
            self.data //= divider

    def __values(self):
        return self.data[0, 0], self.data[0, 1], self.data[1, 0], self.data[1, 1]


class GeneralizedContinuedFraction(object):
    def __init__(self, a_=None, b_=None) -> None:
        """
        each matrix of the series will be:
        M_i =   [[ 0, b_i],
                 [ 1, a_i]]
        this will result in the continued fraction: a_0 +    b_0
                                                          --------
                                                          a_1+ b_1
                                                              -----
                                                                 ..
                                                                  .
        :param a_: a_n series
        :param b_: b_n series
        """
        super().__init__()
        self.mobius = MobiusTransform()
        self.a_ = []
        self.b_ = []
        if (a_ is not None) and (b_ is not None):
            self.extend(a_, b_)

    def evaluate(self, x=None):
        return self.mobius(x) + self.a_[0]

    def extend(self, a_, b_):
        """
        add depth to existing GCF
        """
        self.a_ = (self.a_ + a_).copy()
        self.b_ = (self.b_ + b_).copy()
        for i in range(min(len(a_)-1, len(b_))):
            mat = np.array([[0, b_[i]], [1, a_[i+1]]], dtype=object)
            self.mobius *= MobiusTransform(mat)

    def __len__(self, item):
        return len(self.a_)

    def sym_expression(self, n):
        x = symbols('..')
        eq = x
        for i in reversed(range(n)):
            eq = self.a_[i] + self.b_[i] / eq
        return eq

    def print(self, n=3):
        pprint(self.sym_expression(n))


class SimpleContinuedFraction(GeneralizedContinuedFraction):
    def __init__(self, const_gen, n):
        """
        regular continued fraction
        :param const_gen: must be a generator function (implemented const_gen()). this will give us the constant
        :type const_gen: function
        :param n: depth of CF
        """
        a_ = create_simple_continued_fraction(const_gen, n)
        super().__init__(a_, [1] * len(a_))

    def __str__(self):
        return str(self.a_)


def create_simple_continued_fraction(const_gen, n):
    """
    the known algorithm written with mobius transforms:
        1) tmp = 1/x
        2) a_i = floor(tmp)
        3) x = 1/x - a_i
    instead of calculating on x, we calculate the transforms along the way on x.
    :param const_gen: generator for constant to make continued fraction of
    :param n: depth of continued fraction.
    :return: simple continued fraction of const, with depth n.
    TODO - wasteful to calculate beyond first decimal point.. in (**)
    """
    k = MobiusTransform()
    a_ = [np.floor(const_gen(), dtype=object)]
    const = const_gen() - a_[0]    # could be useful to have better precision along the way
    for i in range(1, n):
        k_rcp = MobiusTransform(np.array([[0, 1], [1, 0]], dtype=object)) * k     # 1) calculate floor(1/x)
        rcp = k_rcp(const)                                          # 1) (**)
        a_.append(floor(rcp))                                       # 2) find a_i
        arr = np.array([[-a_[i], 1], [1, 0]], dtype=object)         # 3) x = 1/x - a_i
        k = MobiusTransform(arr) * k                                # 3) in fact, x = [[1, -a_i], [1, 0]] on
    return a_


def check_and_modify_precision(const, transform, const_gen, offset):
    """
    tried to use this to calculate and enlarge precision along the way.. but it only made it slower.
    right now this is unused.
    :param const:
    :param transform:
    :param const_gen:
    :param offset:
    :return:
    """
    const_work = const
    while True:
        epsilon = dec(2) ** (-mpmath.mp.prec + 5)
        const_d = const_work - epsilon
        const_u = const_work + epsilon
        next_d = floor(transform(const_d))
        next_u = floor(transform(const_u))
        if next_d == next_u:
            return next_d, const_work
        mpmath.mp.prec += 100
        const_work = const_gen() + offset
