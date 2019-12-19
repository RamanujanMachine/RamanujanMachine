from collections import namedtuple
from sympy import pi as pi
from sympy import E as e
from sympy import zeta

MasseySeries = namedtuple('MasseySeries', 'shift_reg initials')
CFData = namedtuple('CFData', 'lhs rhs_an rhs_bn')

pi_cf = {
    'pi1': CFData(4 / pi, MasseySeries([1, -2, 1], [1, 3]), MasseySeries([1, -3, 3, -1], [1, 4, 9])),
    # THESE ARE COMMENTED OUT BECAUSE THE CONVERGE TOO SLOW (sublinear from wikipedia)
    # 'pi2': CFData(lambda: pi+3, MasseySeries([1, -1], [6]),         MasseySeries([1, -3, 3, -1], [1, 9, 25])),
    # 'pi3': CFData(lambda: (4/pi)+1, MasseySeries([1, -1], [2]),     MasseySeries([1, -3, 3, -1], [1, 9, 25])),
    'pi4': CFData(-8 / pi, MasseySeries([1, -2, 1], [-2, -6]), MasseySeries([1, -3, 3, -1], [4, 16, 36])),
    'pi5': CFData(1 + 4 / (pi - 4), MasseySeries([1, -2, 1], [-4, -7]),
                  MasseySeries([1, -3, 3, -1], [-2, -9, -20])),
    'pi6': CFData(4 / (3 * pi - 8), MasseySeries([1, -2, 1], [3, 6]),
                  MasseySeries([1, -3, 3, -1], [-1, -6, -15])),
    'pi7': CFData(2 / (pi + 2), MasseySeries([1, -2, 1], [0, 3]), MasseySeries([1, -3, 3, -1], [1, -2, -9])),
    'pi8': CFData(4 / pi, MasseySeries([1, -2, 1], [1, 4]), MasseySeries([1, -3, 3, -1], [1, -2, -9])),
    'pi9': CFData(4 / (2 - pi), MasseySeries([1, -2, 1], [-3, -5]), MasseySeries([1, -3, 3, -1], [3, 8, 15])),
}

e_cf = {
    'e1': CFData(e / (e - 2), MasseySeries([1, 0, 0, -1, 0, 0, -1, 0, 0, 1], [3, 1, 3, 1, 1, 1, 3, 3, 3]),
                 MasseySeries([1, -1], [1])),
    'e2': CFData(e / (e - 2), MasseySeries([1, -2, 1], [4, 5]), MasseySeries([1, -2, 1], [-1, -2])),
    'e3': CFData(e, MasseySeries([1, -2, 1], [3, 4]), MasseySeries([1, -2, 1], [-1, -2])),
    'e4': CFData(e - 1, MasseySeries([1, 0, 0, -2, 0, 0, 1], [1, 1, 2, 1, 1, 4]), MasseySeries([1, -1], [1])),
}

zeta_cf = {
    # THESE ARE COMMENTED OUT BECAUSE THE CONVERGE TOO SLOW
    # 'zeta3': CFData(lambda: 1/zeta(3), MasseySeries([1, -4, 6, -4, 1], [0**3+1**3, 1**3+2**3, 2**3+3**3, 3**3+4**3]),
    #             MasseySeries([1, -7, 21, -35, 35, -21, 7, -1], [-1**6, -2**6, -3**6, -4**6, -5**6, -6**6, -7**6]))
    #
    # 'zeta4': CFData(lambda: 13 / zeta(4),
    #                # https://tpiezas.wordpress.com/2012/05/06/zudilins-continued-fraction-for-zeta4/
    #                MasseySeries([1, -6, 15, -20, 15, -6, 1], [12, 2142, 26790, 142968, 500688, 1363362]),
    #                MasseySeries([1, -9, 36, -84, 126, -126, 84, -36, 9, -1],
    #                             [24, 26880, 1574640, 28114944, 262500000, 1627547904, 7609537320, 28940697600,
    #                              94014038664])),
    'zeta3_2': CFData(5 / (2 * zeta(3)),
                      MasseySeries([1, -4, 6, -4, 1], [2 + 0 * 2 * 4, 2 + 1 * 3 * 7, 2 + 2 * 4 * 10, 2 + 3 * 5 * 13]),
                      MasseySeries([1, -7, 21, -35, 35, -21, 7, -1],
                                   [2 * (1 ** 5) * 1, 2 * (2 ** 5) * 3, 2 * (3 ** 5) * 5, 2 * (4 ** 5) * 7,
                                    2 * (5 ** 5) * 9, 2 * (6 ** 5) * 11, 2 * (7 ** 5) * 13])),

    # https://tpiezas.wordpress.com/2012/05/04/continued-fractions-for-zeta2-and-zeta3/
    'zeta3_3': CFData(6 / zeta(3),
                      MasseySeries([1, -4, 6, -4, 1], [5, 117, 535, 1463]),
                      MasseySeries([1, -7, 21, -35, 35, -21, 7, -1],
                                   [-1**6, -2**6, -3**6, -4**6, -5**6, -6**6, -7**6])),
    'zeta2': CFData(5 / zeta(2),
                    MasseySeries([1, -3, 3, -1],
                                 [3 - 11 * 1 + 11 * (1 ** 2), 3 - 11 * 2 + 11 * (2 ** 2), 3 - 11 * 3 + 11 * (3 ** 2)]),
                    MasseySeries([1, -5, 10, -10, 5, -1], [1 ** 4, 2 ** 4, 3 ** 4, 4 ** 4, 5 ** 4])),
}

weird_stuff = {
    'one': CFData(-1.0, MasseySeries([1, -2, 1], [-4, -7]), MasseySeries([1, -3, 3, -1], [-9, -20, -35])),
}

# TODO: Q1) are all LSFR binomial coefficients with alternating signs (X-Y)^d?
#       Answer - the solution of these recurrence equations are polynomials (of degree d-1).
#       since the ramanujan project looked for polynomials, of course the recurrence equations are of this form.
# TODO: Q2) maybe some connection between the generating function of the recurrence equation to the the LHS?
#          (probably not...)
# TODO: Q3) all of pi formulas are of the same LSFR?
# TODO: Q4) does every representation of e fit in a RFC?
# TODO: Q5) are all rhs are mobius transforms on the constant? maybe implicit graph search will find us the correct one?
# TODO: 6) write function to determine convergence rate.
