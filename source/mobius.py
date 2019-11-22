import numpy as np
from math import gcd
from mpmath import mpf as dec


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

    def __call__(self, x):
        """
        apply transformation on X.
        :param x (for our project it will always be 1):
        :return: result of mobius transformation
        """
        a, b, c, d = self.__values()
        numerator = dec(a*x + b)
        denominator = dec(c*x + d)
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


class ContinuedFraction(object):
    def __init__(self, a_, b_) -> None:
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
        self.a0 = a_[0]
        for i in range(min(len(a_)-1, len(b_))):
            mat = np.array([[0, b_[i]], [1, a_[i+1]]])
            self.mobius *= MobiusTransform(mat)

    def evaluate(self):
        return self.mobius(1) + self.a0
