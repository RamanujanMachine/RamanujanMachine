from massey import create_series_from_shift_reg
from collections import namedtuple
from mpmath import pi as pi
from mpmath import e as e
from mpmath import zeta

MasseySeries = namedtuple('MasseySeries', 'shift_reg initials')
CFData = namedtuple('CFData', 'lhs rhs_an rhs_bn')

pi_cf = {
    'pi1': CFData(lambda: 4 / pi, MasseySeries([1, -2, 1], [1, 3]), MasseySeries([1, -3, 3, -1], [1, 4, 9])),
    # THESE ARE COMMENTED OUT BECAUSE THE CONVERGE TOO SLOW (sublinear from wikipedia)
    # 'pi2': CFData(lambda: pi+3, MasseySeries([1, -1], [6]),         MasseySeries([1, -3, 3, -1], [1, 9, 25])),
    # 'pi3': CFData(lambda: (4/pi)+1, MasseySeries([1, -1], [2]),     MasseySeries([1, -3, 3, -1], [1, 9, 25])),
    'pi4': CFData(lambda: -8 / pi, MasseySeries([1, -2, 1], [-2, -6]), MasseySeries([1, -3, 3, -1], [4, 16, 36])),
    'pi5': CFData(lambda: 1 + 4 / (pi - 4), MasseySeries([1, -2, 1], [-4, -7]),
                  MasseySeries([1, -3, 3, -1], [-2, -9, -20])),
    'pi6': CFData(lambda: 4 / (3 * pi - 8), MasseySeries([1, -2, 1], [3, 6]),
                  MasseySeries([1, -3, 3, -1], [-1, -6, -15])),
    'pi7': CFData(lambda: 2 / (pi + 2), MasseySeries([1, -2, 1], [0, 3]), MasseySeries([1, -3, 3, -1], [1, -2, -9])),
    'pi8': CFData(lambda: 4 / pi, MasseySeries([1, -2, 1], [1, 4]), MasseySeries([1, -3, 3, -1], [1, -2, -9])),
    'pi9': CFData(lambda: 4 / (2 - pi), MasseySeries([1, -2, 1], [-3, -5]), MasseySeries([1, -3, 3, -1], [3, 8, 15])),
}

e_cf = {
    'e1': CFData(lambda: e / (e - 2), MasseySeries([1, 0, 0, -1, 0, 0, -1, 0, 0, 1], [3, 1, 3, 1, 1, 1, 3, 3, 3]),
                 MasseySeries([1, -1], [1])),
    'e2': CFData(lambda: e / (e - 2), MasseySeries([1, -2, 1], [4, 5]), MasseySeries([1, -2, 1], [-1, -2])),
    'e3': CFData(lambda: e, MasseySeries([1, -2, 1], [3, 4]), MasseySeries([1, -2, 1], [-1, -2])),
    'e4': CFData(lambda: e - 1, MasseySeries([1, 0, 0, -2, 0, 0, 1], [1, 1, 2, 1, 1, 4]), MasseySeries([1, -1], [1])),
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
    'zeta3_2': CFData(lambda: 5 / (2 * zeta(3)),
                      MasseySeries([1, -4, 6, -4, 1], [2 + 0 * 2 * 4, 2 + 1 * 3 * 7, 2 + 2 * 4 * 10, 2 + 3 * 5 * 13]),
                      MasseySeries([1, -7, 21, -35, 35, -21, 7, -1],
                                   [2 * (1 ** 5) * 1, 2 * (2 ** 5) * 3, 2 * (3 ** 5) * 5, 2 * (4 ** 5) * 7,
                                    2 * (5 ** 5) * 9, 2 * (6 ** 5) * 11, 2 * (7 ** 5) * 13])),

    # https://tpiezas.wordpress.com/2012/05/04/continued-fractions-for-zeta2-and-zeta3/
    'zeta3_3': CFData(lambda: 6 / zeta(3),
                      MasseySeries([1, -4, 6, -4, 1], [5, 117, 535, 1463]),
                      MasseySeries([1, -7, 21, -35, 35, -21, 7, -1],
                                   [-1**6, -2**6, -3**6, -4**6, -5**6, -6**6, -7**6])),
    'zeta2': CFData(lambda: 5 / zeta(2),
                    MasseySeries([1, -3, 3, -1],
                                 [3 - 11 * 1 + 11 * (1 ** 2), 3 - 11 * 2 + 11 * (2 ** 2), 3 - 11 * 3 + 11 * (3 ** 2)]),
                    MasseySeries([1, -5, 10, -10, 5, -1], [1 ** 4, 2 ** 4, 3 ** 4, 4 ** 4, 5 ** 4])),
}


# TODO: 1) are all LSFR binom coefficients with alternating signs?
# TODO: 2) notice all shift register are symmetric!
# TODO: 3) shift register -> recursion -> generating functions -> ???
# TODO: 4) all of pi are the same?
# TODO: 5) does every representation of e fit in a RFC?
