import math
import mpmath

from .RelativeGCFEnumerator import RelativeGCFEnumerator
from collections import namedtuple

CONVERGENCE_THRESHOLD = 0.1
BURST_NUMBER = 200
FIRST_ENUMERATION_MAX_DEPT = 5_000
MIN_ITERS = 1

Match = namedtuple('Match', 'rhs_an_poly rhs_bn_poly')
RefinedMatch = namedtuple('RefinedMatch', 'rhs_an_poly rhs_bn_poly val c_top c_bot precision')


def check_for_fr(an_iterator, bn_iterator, an_deg, burst_number=BURST_NUMBER, min_iters=MIN_ITERS):
    """
    As the calculation for p and q goes on, the GCD for the two grows. 
    We've noticed that conjectures tends to have a GCD that grows in a super exponential manner (we call that Factorial
    Reduction).
    This function test if a GCF has factorial reduction.
    """
    calculated_values = []
    num_of_calculated_vals = 0

    prev_q = 0
    q = 1
    prev_p = 1
    # This is a ugly hack but it works. a[0] is handled before the rest here:
    p = an_iterator.__next__()  # will place a[0] to p
    bn_iterator.__next__()  # b0 is discarded

    next_gcd_calculation = burst_number if burst_number >= min_iters else min_iters

    for i, (a_i, b_i) in enumerate(zip(an_iterator, bn_iterator)):
        tmp_a = q
        tmp_b = p

        q = a_i * q + b_i * prev_q
        p = a_i * p + b_i * prev_p

        prev_q = tmp_a
        prev_p = tmp_b

        if i == next_gcd_calculation:
            num_of_calculated_vals += 1
            next_gcd_calculation += burst_number

            calculated_values.append(
                mpmath.log(mpmath.mpf(math.gcd(p, q))) / mpmath.mpf(i) +
                an_deg * (-mpmath.log(i) + 1)
            )

            # This test fails for GCF without FR. Checking it early allows us do discard a lot of GCFs
            if num_of_calculated_vals >= 3 and \
                    abs(calculated_values[-2] - calculated_values[-1]) > \
                    abs(calculated_values[-2] - calculated_values[-3]):
                return False, i

            if num_of_calculated_vals >= 2 and \
                    abs(calculated_values[-2] - calculated_values[-1]) < CONVERGENCE_THRESHOLD:
                return True, i

    return False, i


class FREnumerator(RelativeGCFEnumerator):
    """
    This enumerator checks the Factorial Reduction property of GCFs as the first step of the enumeration.
    In the FR test we don't compute the GCF's value, or even compare it to a LHS.
    Fractions that have FR will be computed to a higher dept using RelativeGCFEnumerator's implementation.
    The computed values are then fed into a PSLQ that tries to find a suitable LHS.
    """

    def __init__(self, *args, **kwargs):
        print('checking for FR enumerator')
        super().__init__(None, *args, **kwargs)

    def _first_enumeration(self, print_results: bool):
        """
        Test all GCFs in the domain for FR.
        """
        results = []  # list of intermediate results        
        all_items_calculated = []
        for an_iter, bn_iter, metadata in self._iter_domains_with_cache(FIRST_ENUMERATION_MAX_DEPT):
            has_fr, items_calculated = check_for_fr(an_iter, bn_iter, self.poly_domains.get_an_degree(metadata.an_coef))
            if has_fr:
                all_items_calculated.append(items_calculated)
                if print_results:
                    print(f"found a GCF with FR:\n\tan: {metadata.an_coef}\n\tbn: {metadata.bn_coef}")
                # Key is useless here :)
                results.append(Match(metadata.an_coef, metadata.bn_coef))

        return results

    def _improve_results_precision(self, intermediate_results, verbose=True):
        """
        Calculates GCFs to a higher dept using RelativeGCFEnumerator's implementation.
        We then feed those results and the constant given to a PSLQ, that tries to find a suitable LHS.

        Notice-
        The second part of this function (PSLQ), logically belongs to the next step of the algorithm - the 
        result refinement part. It is implemented here, because the next function is not parallelized over
        different processes or clients, and we want the PSLQ to be parallelized as well. 
        """
        precise_intermediate_results = super()._improve_results_precision(intermediate_results, verbose)

        pslq_results = []
        consts = [i() for i in self.constants_generator]
        for match, val, precision in precise_intermediate_results:
            mpf_val = mpmath.mpf(val)
            pslq_res = mpmath.pslq(
                [1, consts[0], consts[0] * consts[0], -mpf_val, -consts[0] * mpf_val, -consts[0] * consts[0] * mpf_val],
                tol=10 ** (1 - precision))
            if pslq_res:
                pslq_results.append(RefinedMatch(*match, val, pslq_res[:3], pslq_res[3:], precision))
            else:
                pslq_results.append(RefinedMatch(*match, val, None, None, precision))

        return pslq_results

    def _refine_results(self, intermediate_results, print_results=True):
        return intermediate_results
