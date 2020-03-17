import numpy as np
from math import gcd, floor, ceil
from mpmath import mpf as dec
import mpmath
from sympy import Symbol, pprint
# from ortools.linear_solver.pywraplp import Solver


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

    def sym_expression(self, x):
        """
        get the symbolic expression of the transformation.
        :param x: expression to use as the operand of the transformation
        """
        a, b, c, d = self.__values()
        return (a*x + b) / (c*x + d)

    def pprint(self, x=Symbol('x')):
        """
        pretty print the mobius transform.
        :param x: (optional) expression to print as the operand of the transformation
        """
        sym = self.sym_expression(x)
        pprint(sym)

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
            numerator = dec(a * x + b)
            denominator = dec(c * x + d)
        else:
            numerator = dec(b)
            denominator = dec(d)
        return numerator / denominator

    def __eq__(self, other):
        """
        Compare with another mobius.
        :param other: another mobius object.
          :return: True if equal. Else False.
        """
        if not isinstance(other, MobiusTransform):
            raise TypeError("Comparision of wrong types")
        if np.array_equal(self.data, other.data):
            return True
        else:
            return False

    def normalize(self):
        """
        to be called after every operation. this should prevent coefficients from exploding.
        """
        a, b, c, d = self.__values()
        divider = gcd(gcd(a, b), gcd(c, d))
        if divider != 1:
            self.data //= divider

    def reciprocal(self):
        """
        let T be the transform, and x be the operand. the output transform does 1/T(x).
        :return: a reciprocal transform
        """
        a, b, c, d = self.__values()
        tmp = MobiusTransform(np.array([[c, d], [a, b]], dtype=object))
        return tmp

    def inverse(self):
        """
        :return: the inverse transformation
        """
        a, b, c, d = self.__values()
        det = a * d - b * c
        tmp = MobiusTransform(det * np.array([[d, -b], [-c, a]], dtype=object))
        tmp.normalize()
        return tmp

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

    def evaluate(self):
        """
        evaluate numerically the convergent of the the GCF
        :return: mpmath decimal.
        """
        return self.mobius(None) + self.a_[0]

    def extend(self, a_, b_):
        """
        add depth to existing GCF
        """
        self.a_ = (self.a_ + a_).copy()
        self.b_ = (self.b_ + b_).copy()
        for i in range(min(len(a_) - 1, len(b_))):
            mat = np.array([[0, b_[i]], [1, a_[i + 1]]], dtype=object)
            self.mobius *= MobiusTransform(mat)

    def __len__(self, item):
        return len(self.a_)

    def sym_expression(self, n):
        """
        convert convergent no. n into a sympy expression (non simplified)
        this is used to pretty print the GCF
        :param n: depth of convergent
        :return: sym expression
        """
        x = Symbol('..')
        eq = x
        for i in reversed(range(n)):
            eq = self.a_[i] + self.b_[i] / eq
        return eq

    def print(self, n=3):
        """
        pretty print the GCF
        :param n: depth to print
        """
        pprint(self.sym_expression(n))

    @classmethod
    def from_irrational_constant(cls, const_gen, b_):
        """
        if b_ is all '1's and s_ is not defined, this is the known algorithm written with mobius transforms:
            1) tmp = 1/x
            2) a_i = floor(tmp)
            3) x = 1/x - a_i
        instead of calculating on x, we calculate the transforms along the way on x.
        TODO - wasteful to calculate beyond first decimal point.. in (**)
        :param const_gen: must be a generator function (implemented const_gen()). this will give us the constant
        :type const_gen: function
        :param b_: series of nominators for the generalized continued fraction
        """
        const = const_gen()  # could be useful to have better precision along the way
        a_ = [floor(const) if b_[0] > 0 else ceil(const)]
        k = MobiusTransform(np.array([[1, -a_[0]], [0, 1]]))  # x = x - a[0]
        for i in range(1, len(b_)):
            k_rcp = MobiusTransform(np.array([[0, b_[i - 1]], [1, 0]], dtype=object)) * k  # 1) calculate floor(b[i]/x)
            try:
                rcp = k_rcp(const)  # 1) (**)
            except ZeroDivisionError:
                print("create simple continued fraction finished sooner than expected resulting in finite fraction\n"
                      "this may be due to a rational number given as input, or insufficient precision")
                # return cls(a_, b_)
                raise ZeroDivisionError
            a_.append(floor(rcp) if b_[i] > 0 else ceil(rcp))  # 2) find a_i
            next_transform = MobiusTransform(np.array([[0, b_[i-1]], [1, a_[i]]], dtype=object))  # 3) x = b[i]/x - a[i]
            k = next_transform.inverse() * k
        return cls(a_, b_)

    def __eq__(self, other):
        if not isinstance(other, GeneralizedContinuedFraction):
            raise TypeError("Comparision of wrong types")
        for i in range(min(len(self.a_), len(other.a_))):
            res = (self.a_[i] == other.a_[i]) and (self.b_[i] == other.b_[i])
            if not res:
                return False
        return True


class SimpleContinuedFraction(GeneralizedContinuedFraction):
    def __init__(self, a_=None):
        """
        regular continued fraction
        :param a_: must be a generator function (implemented const_gen()). this will give us the constant
        """
        if a_ is None:
            super().__init__()
        else:
            super().__init__(a_, [1] * len(a_))

    def __str__(self):
        return str(self.a_)

    @classmethod
    def from_irrational_constant(cls, const_gen, n, **kwargs):
        """
        the known algorithm written with mobius transforms:
            1) tmp = 1/x
            2) a_i = floor(tmp)
            3) x = 1/x - a_i
        instead of calculating on x, we calculate the transforms along the way on x.
        TODO - wasteful to calculate beyond first decimal point.. in (**)
        :param const_gen: must be a generator function (implemented const_gen()). this will give us the constant
        :type const_gen: function
        :param n: depth of CF
        """
        gcf = GeneralizedContinuedFraction.from_irrational_constant(const_gen, [1] * n)
        return cls(gcf.a_)


class EfficientGCF(object):
    def __init__(self, a_, b_):
        """
        efficient calculation of general continued fraction - only build and evaluate.
        the calculation is done with the recursive formula of the gcf convergent.
        ref: https://en.wikipedia.org/wiki/Generalized_continued_fraction (switch names between a and b in ref)
        :param a_: an series
        :param b_: bn series
        """
        self.prev_A = 0
        self.A = 1
        self.prev_B = 1
        self.B = a_[0]

        for i in range(1, len(a_)):
            tmp_a = self.A
            tmp_b = self.B
            self.A = a_[i] * self.A + b_[i - 1] * self.prev_A
            self.B = a_[i] * self.B + b_[i - 1] * self.prev_B
            self.prev_A = tmp_a
            self.prev_B = tmp_b

    def evaluate(self):
        if self.A == 0:
            return dec(0)
        return dec(self.B) / dec(self.A)


# def find_transform(x, y, limit, threshold=1e-7):
#     """
#     find a integer solution to ax + b - cxy - dy = 0
#     this will give us the mobius transform: T(x) = y
#     :param x: numeric constant to check
#     :param y: numeric manipulation of constant
#     :param limit: range to look at
#     :param threshold: optimal solution threshold.
#     :return MobiusTransform in case of success or None.
#     """
#     x1 = x
#     x2 = dec(1.0)
#     x3 = -x*y
#     x4 = -y
#     solver = Solver('mobius', Solver.CBC_MIXED_INTEGER_PROGRAMMING)
#     a = solver.IntVar(-limit, limit, 'a')
#     b = solver.IntVar(-limit, limit, 'b')
#     c = solver.IntVar(-limit, limit, 'c')
#     d = solver.IntVar(-limit, limit, 'd')
#     f = solver.NumVar(0, 1, 'f')
#     solver.Add(f == (a*x1 + b*x2 + c*x3 + d*x4))
#     solver.Add(a*x1 + b >= 1)   # don't except trivial solutions and remove some redundancy
#     solver.Minimize(f)
#     status = solver.Solve()
#     if status == Solver.OPTIMAL:
#         if abs(solver.Objective().Value()) <= threshold:
#             res_a, res_b, res_c, res_d = int(a.solution_value()), int(b.solution_value()),\
#                                          int(c.solution_value()), int(d.solution_value())
#             ret = MobiusTransform(np.array([[res_a, res_b], [res_c, res_d]], dtype=object))
#             ret.normalize()
#             return ret
#     else:
#         return None


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
