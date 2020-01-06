import numpy as np
import math
import mobius
from massey import create_series_from_polynomial, create_series_from_shift_reg, create_series_from_compact_poly
def max_depth_for_n_bits(n):
    depth = 10
    over_n = False
    while (over_n == False):
        rhs_an = create_series_from_polynomial([3, 3], depth)
        rhs_bn = create_series_from_shift_reg([1, -3, 3, -1], [-1 * 1, -2 * 3, -3 * 5], depth)
        GCF = mobius.GeneralizedContinuedFraction(rhs_an, rhs_bn)
        data = GCF.mobius.data
        vals = [int(np.abs(data[0][0])), int(np.abs(data[0][1])), int(np.abs(data[1][0])), int(np.abs(data[1][1]))]
        Max = np.max(vals)
        bits = math.log2(Max)
        if bits > n:
            over_n = True
        depth += 1
    return depth




if __name__ == "__main__":
    depth = max_depth_for_n_bits(64)
    print("Max depth without exceeding {} bits per matrix element is:  {}".format(64, depth))
