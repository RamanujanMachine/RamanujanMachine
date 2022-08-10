import matplotlib.pyplot as plt
import mpmath
from mpmath import mpf as dec
from ramanujan.utils.mobius import GeneralizedContinuedFraction


def calculate_convergence(
    gcf: GeneralizedContinuedFraction, reference, plot=False, title=""
):
    """
    calculate convergence rate of General Continued Fraction (in reference to some constant x).
    the result is the average number of decimal digits, per term of the general continued fraction.
    :param plot: whether or not to plot the graph of the log10 difference between convergent and x.
    :param gcf: General Continued Fraction to calculate.
    :param title: (optional) title of graph.
    :param reference: x
    """
    q_ = [0, 1]
    p_ = [1, gcf.a_[0]]
    log_diff = [mpmath.log10(abs((dec(p_[1]) / dec(q_[1])) - reference))]
    length = min(200, len(gcf.b_))
    for i in range(1, length):
        q_.append(gcf.a_[i] * q_[i] + gcf.b_[i - 1] * q_[i - 1])
        p_.append(gcf.a_[i] * p_[i] + gcf.b_[i - 1] * p_[i - 1])
        if q_[i + 1] == 0:
            part_convergent = 0
        else:
            part_convergent = dec(p_[i + 1]) / dec(q_[i + 1])
        if not mpmath.isfinite(part_convergent):
            length = i
            break
        log_diff.append(mpmath.log10(abs(part_convergent - reference)))
    if plot:
        plt.plot(range(length), log_diff)
        plt.title(title)
        plt.show()
    log_slope = 2 * (log_diff[length - 1] - log_diff[length // 2]) / length
    return -log_slope
