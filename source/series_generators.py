from itertools import chain, product, combinations_with_replacement
from typing import List, Iterator, Callable
from abc import ABC, abstractmethod
from functools import partial
from math import factorial


def number_of_cartesian_product_elements(permutation_options: List[List]):
    """
    returns the number of permutations after calling product(*permutation_options)
    """
    res = 1
    for l in permutation_options:
        res *= len(l)
    return res


def binomial(n, k):
    """
    return ( n )  binomial coefficient
           ( k )
    """
    try:
        binom = factorial(n) // factorial(k) // factorial(n - k)
    except ValueError:
        binom = 0
    return binom


class SeriesGeneratorClass(ABC):
    """ Abstract class for series generator class """

    @abstractmethod
    def get_function(self):
        """
        this will be used to get the series generating function.
        this function must return series generating function of type f(x_, n)
        where x_ is the current permutation, and n is the length of the list
        """
        pass

    @abstractmethod
    def get_num_iterations(self, poly_options):
        """
        return the theoretical number of permutations being iterated
        :param poly_options: coefficient options.
        """
        pass

    @abstractmethod
    def get_iterator(self, poly_options):
        """
        return the permutations iterator
        :param poly_options: coefficient options.
        """
        pass

    help_string = 'this is a help string to be presented in API'


class CartesianProductAnGenerator(SeriesGeneratorClass):
    """
    Default class of polynomial Cartesian product series generator.
    """

    def __init__(self):
        self.function = create_series_from_compact_poly

    def get_function(self) -> Callable[[List[int], int], List[int]]:
        return self.function

    def get_num_iterations(self, poly_a: List[List[int]]):
        #return 2 * number_of_cartesian_product_elements(poly_a)

        # TODO - acutally make the following be:
        # allowing all elements to be negetive beside the leading coef (to save on duplicate runs)
        return number_of_cartesian_product_elements(poly_a[1:]) * len(poly_a[1])

    def get_iterator(self, poly_a: List[List[int]]) -> Iterator:
        """
        return cartesian product iterator
        if poly_a = [[1,2,3],[5,6]], permutations will be { [1,5] ; [1,6] ; [2,5] ; [2,6] ; [3,5] ; [3,6] }
        """
        neg_poly_a = [poly_a[0], *[[-i for i in a] for a in poly_a[1:]]]  # for b_n include negative terms
        return chain(product(*poly_a), product(*neg_poly_a))
        #return product(*poly_a)

    help_string = 'a[n] = n(n(...(x[1]*n + x[0]) + x[2]) + ...) + x[k]. this is the default generator'


class CartesianProductBnGenerator(SeriesGeneratorClass):
    help_string = 'b[n] = n(n(...(x[1]*n + x[0]) + x[2]) + ...) + x[k]. this is the default generator'

    def __init__(self):
        self.function = create_series_from_compact_poly

    def get_function(self) -> Callable[[List[int], int], List[int]]:
        return self.function

    def get_num_iterations(self, poly_b: List[List[int]]):
        return 2 * number_of_cartesian_product_elements(poly_b)

    def get_iterator(self, poly_b: List[List[int]]) -> Iterator:
        """
        return cartesian product iterator including negative options
        if poly_b = [[1,2,3],[5,6]], permutations will be { [1,5] ; [1,6] ; [2,5] ; [2,6] ; [3,5] ; [3,6] }
                                                      +   { [-1,-5] ; [-1,-6] ; [-2,-5] ; [-2,-6] ; [-3,-5] ; [-3,-6] }
        """
        neg_poly_b = [[-i for i in b] for b in poly_b]  # for b_n include negative terms
        #all_possibilites = neg_poly_b + poly_b
        all_possibilites = [[-i for i in b][1:] + b for b in poly_b]

        return product(*all_possibilites)



class CartesianProductBnShift1(CartesianProductBnGenerator):
    help_string = 'b[n] = m(m(...(x[1]*m + x[0]) + x[2]) + ...) + x[k], where m=n+1'

    def __init__(self):
        super().__init__()
        self.function = create_series_from_compact_poly_with_shift1


class CartesianProductAnShift1(CartesianProductAnGenerator):
    help_string = 'a[n] = m(m(...(x[1]*m + x[0]) + x[2]) + ...) + x[k], where m=n+1'

    def __init__(self):
        super().__init__()
        self.function = create_series_from_compact_poly_with_shift1


class CartesianProductAnOddOnly(CartesianProductAnGenerator):
    help_string = 'same as default generator with odd coefficients only'

    def __init__(self):
        super().__init__()

    def get_iterator(self, poly_a: List[List[int]]) -> Iterator:
        odd_poly = [[i for i in a if (i % 2) == 1] for a in poly_a]  # for b_n include negative terms
        return product(*odd_poly)

    def get_num_iterations(self, poly_a: List[List[int]]):
        odd_poly = [[i for i in a if (i % 2) == 1] for a in poly_a]
        return super().get_num_iterations(odd_poly)


class CartesianProductBnShift2n1(CartesianProductBnGenerator):
    help_string = 'b[n] = m(m(...(x[1]*m + x[0]) + x[2]) + ...) + x[k], where m=2*n+1'

    def __init__(self):
        super().__init__()
        self.function = create_series_from_compact_poly_with_shift2n1


class CartesianProductBnCatalan(CartesianProductBnGenerator):
    help_string = 'x[0]*(2*n+1)^4 + x[1]*(2*n+1)^3'

    def __init__(self):
        super().__init__()
        self.function = catalan_bn_generator


class CartesianProductZetaBn(CartesianProductBnGenerator):
    help_string = 'b[n] = x[0]*(n+1)^d - x[1]*(n+1)^(d-1).  where d=function_value\n' \
                  'this was found to be useful for zeta values searches.'

    def __init__(self, zeta_value):
        super().__init__()
        self.function = partial(create_zeta_bn_series, zeta_value * 2)


class CartesianProductZeta3An(CartesianProductAnGenerator):
    help_string = 'Generator3[x3_, x0_] := {x0, 2 *x0 + x3, 3*x3, 2*x3}.\n' \
                  'this was found to be useful for zeta3 searches'

    def __init__(self):
        super().__init__()
        self.function = zeta3_an_generator

class CartesianProductZeta3N6Bn(CartesianProductBnGenerator):
    help_string = 'Generator3[x0_] := {x_[0]*n^6}.\n' \
                  'this was found to be useful for zeta3 searches'

    def __init__(self):
        super().__init__()
        self.function = zeta3_bn_n6_generator

class CartesianProductZeta3N6ComplementAn(CartesianProductAnGenerator):
    help_string = 'Generator3[x3_, x0_] := {(2*i - 1) * (i * (i - 1) + (x_[0]**2 + 1)/2)}.\n' \
                  'this was found to be useful for zeta3 searches'

    def __init__(self):
        super().__init__()
        self.function = zeta3_an_n6_complement_generator

class CartesianProductZeta5An(CartesianProductAnGenerator):
    help_string = 'Generator5[x5_, x3_, x0_] := {x0, 2 *x0 + x3 - 2 *x5, 3*x3 - 5*x5, 2*x3, 5*x5, 2*x5}'

    def __init__(self):
        super().__init__()
        self.function = zeta5_an_generator


def create_series_from_polynomial(poly_a, n):
    """
    given a polynomial P(x) as generator return An when A[i] = P(i)
    :param poly_a:
    :type poly_a: list
    :param n: number of values to generate
    :type n: int
    :return: series with length of n
    """
    a_ = []
    for i in range(n):
        a_i = 0
        pol = 1
        for j in range(len(poly_a)):
            a_i += pol * poly_a[j]
            pol *= i
        a_.append(a_i)
    return a_


def create_series_from_shift_reg(poly_a, initials, n):
    """
    this is the reversed action to the massey algorithm.
    given a polynomial P as shift register, generate a series.
    series follows the following rule: A[n] = sum i=1 to deg(P): -(P[i] * A[n-i])
    :param poly_a: P(x) : same format as the massey algorithm output (expects poly_a[0] = 1)
    :param initials: starting conditions of the series. must be of length = deg(P)
    :param n: number of values to generate
    :return: series with length of n
    """
    assert ((len(poly_a) - len(initials)) == 1), 'There should be deg(poly) initial conditions'
    assert (len(poly_a) > 0) and (poly_a[0] == 1), 'Illegal polynomial - first coefficient must be 1'
    a_ = []
    for i in range(len(initials)):
        a_.append(initials[i])
    for i in range(len(a_), n):
        a_i = 0
        for j in range(1, len(poly_a)):
            a_i -= poly_a[j] * a_[i - j]
        a_.append(a_i)
    return a_


def create_series_from_compact_poly(poly_a, n):
    """
    create a series of type n(n(...(a[0]*n + a[1]) + a[2]) + ...) + a[k]
    :param poly_a: a[k] coefficients
    :param n: length of series
    :return: a list of numbers in series
    """
    ret = []
    for i in range(n):
        tmp = 0
        for c in poly_a:
            tmp *= i
            tmp += c
        ret.append(tmp)
    return ret

def iter_series_items_from_compact_poly(poly_coef, max_runs, start_n=1):
    """
    create a series of type n(n(...(a[0]*n + a[1]) + a[2]) + ...) + a[k]
    :param poly_a: a[k] coefficients
    :param n: length of series
    :return: a list of numbers in series
    """
    for i in range(start_n, max_runs):
        tmp = 0
        for c in poly_coef:
            tmp *= i
            tmp += c
        yield tmp


def create_series_from_compact_poly_with_shift1(poly_a, n):
    """
    create a series of type m(m(...(a[1]*m + a[0]) + a[2]) + ...) + a[k],
    where m=n+1
    :param poly_a: a[k] coefficients
    :param n: length of series
    :return: a list of numbers in series
    """
    ret = []
    for m in range(1, n + 1):
        tmp = 0
        for c in poly_a:
            tmp *= m
            tmp += c
        ret.append(tmp)
    return ret


def create_series_from_compact_poly_with_shift2n1(poly_a, n):
    """
    create a series of type m(m(...(a[1]*m + a[0]) + a[2]) + ...) + a[k],
    where m=2n+1
    :param poly_a: a[k] coefficients
    :param n: length of series
    :return: a list of numbers in series
    """
    ret = []
    for m in range(1, 2 * n + 1, 2):
        tmp = 0
        for c in poly_a:
            tmp *= m
            tmp += c
        ret.append(tmp)
    return ret


def create_zeta_bn_series(deg, poly_a, n):
    """
    create a series of type: a[0]*(n+1)^d - a[1]*(n+1)^(d-1)
    apparently, this is useful for building zeta function CFs
    :param deg: d
    :param poly_a: coefficients a[0], a[1]. if more exist, they are ignored.
    :param n: depth of series
    :return: a list of numbers in series
    """
    ret = []
    for i in range(n):
        res = poly_a[0] * ((i + 1) ** deg) + poly_a[1] * ((i + 1) ** (deg - 1))
        ret.append(res)
    return ret


def zeta3_an_generator(x_, n):
    """
    Generator3[x3_, x0_] := {x0, 2 *x0 + x3, 3*x3, 2*x3}
    :param x_:
    :param n:
    """
    ret = []
    for i in range(n):
        res = x_[0] + (2 * x_[0] + x_[1]) * i + 3 * x_[1] * (i ** 2) + 2 * x_[1] * (i ** 3)
        ret.append(res)
    return ret


def zeta3_bn_n6_generator(x_, n):
    """
    Generator3_n6[x3_, x0_] := {x0, 2 *x0 + x3, 3*x3, 2*x3}
    :param x_:
    :param n:
    """
    ret = []
    for i in range(1, n):
        res = x_[0] * (i ** 6) 
        ret.append(res)
    return ret


def zeta3_an_n6_complement_generator(x_, n):
    """
    Generator3_n6[x3_, x0_] := {x0, 2 *x0 + x3, 3*x3, 2*x3}
    :param x_:
    :param n:
    """
    print("="*20)
    ret = []
    for i in range(n):
        res = (2*i - 1) * (i * (i - 1) + (x_[0]**2 + 1)/2)
        ret.append(res)
    return ret


def zeta5_an_generator(x_, n):
    """
    Generator5[x5_, x3_, x0_] := {x0, 2 *x0 + x3 - 2 *x5, 3*x3 - 5*x5, 2*x3, 5*x5, 2*x5}
    :param x_:
    :param n:
    """
    ret = []
    for i in range(n):
        res = x_[0] + (2 * x_[0] + x_[1] - 2 * x_[2]) * i + (3 * x_[1] - 5 * x_[2]) * (i ** 2) + \
              2 * x_[1] * (i ** 3) + 5 * x_[2] * (i ** 4) + 2 * x_[2] * (i ** 5)
        ret.append(res)
    return ret


def catalan_bn_generator(x_, n):
    """
    create a series of type: x[0]*(2*n+1)^4 + x[1]*(2*n+1)^3
    :param x_:
    :param n:
    :return:
    """
    ret = []
    for i in range(n):
        res = x_[0] * (2 * i + 1) ** 4 + x_[1] * (2 * i + 1) ** 3
        ret.append(res)
    return ret


def multiples_generator(factors, x_, n):
    """
    options:

    :param factors:
    :param x_:
    :param n:
    """
    ret = []
    x0 = x_[0]
    x = x_[1]
    for i in range(n):
        res = x0
        for k in x:
            res *= factors[k][i]
        ret.append(res)
    return ret


class IntegerFactor(SeriesGeneratorClass):
    def __init__(self, degree, max_length):
        self.deg = degree
        self.factors = [
            [i + 1 for i in range(max_length)],
            [i + 2 for i in range(max_length)],
            [i + 3 for i in range(max_length)],
            [i + 4 for i in range(max_length)],
            [i + 5 for i in range(max_length)],
            [2 * i for i in range(max_length)],
            [2 * i + 1 for i in range(max_length)],
            [2 * i + 3 for i in range(max_length)],
            [2 * i + 5 for i in range(max_length)],
            [2 * i - 1 for i in range(max_length)]
        ]

    def get_num_iterations(self, poly_options):
        poly_real = [i for i in poly_options[0] if i != 0] + [-i for i in poly_options[0] if i != 0]
        return binomial(len(self.factors) + self.deg - 1, self.deg) * len(poly_real)

    def get_iterator(self, poly_options):
        factor_indices = [i for i in range(len(self.factors))]
        permutations = combinations_with_replacement(factor_indices, self.deg)
        poly_real = [i for i in poly_options[0] if i != 0] + [-i for i in poly_options[0] if i != 0]
        return product(poly_real, permutations)

    def get_function(self):
        return partial(multiples_generator, self.factors)

    help_string = 'b[n] = x[0]*(term1)^d1*(term2)^d2*..., where sum(d)=poly_b_order\n' \
                  'this is a unique generator using combination permutations instead of ' \
                  'cartesian product. current terms are: ' + \
                  '(n+1), (n+2), (n+3), (2n), (2n-1), (2n+3), (2n+5)'
