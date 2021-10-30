import math
import mpmath

from .RelativeGCFEnumerator import RelativeGCFEnumerator
from collections import namedtuple
from ramanujan.utils.utils import get_reduced_fraction

CONVERGENCE_THRESHOLD = 0.1
BURST_NUMBER = 200
# This constant is set to ~ 1400, based on results we've found from various executions.
# Changing it may lead to more false positives / missed results.
# We specifically use 1402 to ensure that no items of an/bn are calculated and not used in GCD calculations
FIRST_ENUMERATION_MAX_DEPTH = 1_402
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
    
    p = next(an_iterator)  # will place a[0] to p
    next(bn_iterator)  # b0 is discarded

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

            # The calculated value will converge for GCFs that have FR, but it will not happen monotonicly.
            # We're calculating values once every burst_number iterations, to try and avoid fluctuations' effect
            # If the value still isn't converging to a steady value, we'll halt the calculation early.
            # TODO - add a referrence to Guy & Nadav's paper once its on arxiv
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
    Fractions that have FR will be computed to a higher depth using RelativeGCFEnumerator's implementation.
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
        for an_iter, bn_iter, metadata in self._iter_domains_with_cache(FIRST_ENUMERATION_MAX_DEPTH):
            has_fr, items_calculated = check_for_fr(an_iter, bn_iter, self.poly_domains.get_an_degree(metadata.an_coef))
            if has_fr:
                if print_results:
                    print(f"found a GCF with FR:\n\tan: {metadata.an_coef}\n\tbn: {metadata.bn_coef}")
                # Key is useless here :)
                results.append(Match(metadata.an_coef, metadata.bn_coef))

        return results

    def _improve_results_precision(self, intermediate_results, verbose=True):
        """
        Calculates GCFs to a higher depth using RelativeGCFEnumerator's implementation.
        We then feed those results and the constant given to a PSLQ, that tries to find a suitable LHS.

        Notice-
        The second part of this function (PSLQ), logically belongs to the next step of the algorithm - the 
        result refinement part. It is implemented here, because the next function is not parallelized over
        different processes or clients, and we want the PSLQ to be parallelized as well. 
        """
        precise_intermediate_results = super()._improve_results_precision(intermediate_results, verbose)
        for i in precise_intermediate_results:
            print(i)
        # keeping intermediate values so decedent classes could use this data
        self.precise_intermediate_results = precise_intermediate_results

        print('Running PSLQ')
        pslq_results = []
        # TODO - add docs and change name
        num_items = [1] + [gen() for gen in self.constants_generator]
        num_of_items = len(num_items)

        for match, val, precision in precise_intermediate_results:
            try:
                mpf_val = mpmath.mpf(val)

                denom_items = [-mpf_val * c for c in num_items]
                pslq_res = mpmath.pslq(
                    num_items + denom_items, tol=10 ** (1 - precision),
                    maxcoeff=1_000)

                if pslq_res:
                    # Sometimes, PSLQ can find several results for the same value (e.g. z(3)/(z(3)^2) = 1/z(3))
                    # we'll reduce fraction found to get uniform results
                    reduced_num, reduced_denom = get_reduced_fraction(
                        pslq_res[:num_of_items], pslq_res[num_of_items:], num_of_items - 1)
                    print(f'Found result! an = {match.rhs_an_poly}, bn = {match.rhs_bn_poly}')
                    print(f'Numerator coefs = {reduced_num}, Denominator coefs = {reduced_denom}')
                else:
                    reduced_num, reduced_denom = [], []
            except Exception as e:
                print(f'Exception when using plsq on PCF {match}, {mpmath.nstr(mpf_val, 30)} with constant' +
                      f'{self.const_sym}')
                print(e)
                print('Result saved with None as PSLQ coeffs')
                reduced_num, reduced_denom = None, None

            pslq_results.append(RefinedMatch(*match, val, reduced_num, reduced_denom, precision))

        return pslq_results

    def _refine_results(self, intermediate_results, print_results=True):
        return intermediate_results
