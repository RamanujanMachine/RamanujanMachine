import math
from time import time
import mpmath 

from .RelativeGCFEnumerator import RelativeGCFEnumerator
from .AbstractGCFEnumerator import Match
from ramanujan.constants import g_N_verify_compare_length

CONVERGENCE_THRESHOLD = 0.1
BURST_NUMBER = 200
MIN_ITERS = 1


def check_for_fr(an_iterator, bn_iterator, metadata, an_deg, bn_deg, burst_number=BURST_NUMBER, min_iters=MIN_ITERS):
    calculated_values = []
    num_of_calculated_vals = 0

    prev_q = 0
    q = 1
    prev_p = 1
    # This is a ugly hack but it works. a[0] is handled before the rest here:
    p = an_iterator.__next__() # will place a[0] to p
    bn_iterator.__next__() # b0 is discarded

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
                mpmath.log(mpmath.mpf(math.gcd(p, q))) / mpmath.mpf(i) + \
                    (an_deg) * (-mpmath.log(i) + 1)
            )

            if num_of_calculated_vals >= 3 and \
                abs(calculated_values[-2] - calculated_values[-1]) > abs(calculated_values[-2] - calculated_values[-3]):
                return False, i

            if num_of_calculated_vals >= 2 and \
                abs(calculated_values[-2] - calculated_values[-1]) < CONVERGENCE_THRESHOLD:
                # import ipdb
                # ipdb.set_trace()
                return True , i

    return False, i


class FREnumerator(RelativeGCFEnumerator):
    """
        This enumerator calculates the gcf to an arbitrary dept, until getting to
        a stable value withing the precision required bounds.
        Useful for GCF that converges slowly.
    """

    def __init__(self, *args, **kwargs):
        print('checking for FR enumerator')
        super().__init__(*args, **kwargs)


    def _first_enumeration(self, print_results: bool):
        num_iterations = self.poly_domains.num_iterations

        start = time()
        key_factor = 1 / self.threshold

        results = []  # list of intermediate results        
        all_items_calculated = []
        last_lead = -1001 # random number that will not pop
        for an_iter, bn_iter, metadata in self._iter_domains_with_cache(5000):
            if metadata.an_coef[0] != last_lead:
                last_lead = metadata.an_coef[0]
                print(f"last_lead {last_lead} - {time() - start}")
            
            has_fr, items_calculated = check_for_fr(an_iter, bn_iter, metadata,
                self.poly_domains.get_an_degree(metadata.an_coef),
                self.poly_domains.get_bn_degree(metadata.bn_coef))
            if has_fr:
                all_items_calculated.append(items_calculated)
                print(metadata)
                # Key is useless here :)
                results.append(Match(0, metadata.an_coef, metadata.bn_coef))

        if len(all_items_calculated) != 0:
            print(sum(all_items_calculated)/len(all_items_calculated))

        return results

    def _improve_results_precision(self, intermediate_results, verbose=True):
        precise_intermediate_results = super()._improve_results_precision(intermediate_results, verbose)

        pslq_results = []
        consts = [i() for i in self.constants_generator]
        for match, val, precision in precise_intermediate_results:
            # TODO - create functions that controls number of factors
            mpf_val = mpmath.mpf(val)
            pslq_res = mpmath.pslq(
                [1, consts[0], consts[0]*consts[0], mpf_val, consts[0]*mpf_val, consts[0]*consts[0]*mpf_val], 
                tol=10**(1-precision))
            pslq_results.append((match,val,pslq_res, precision))
            print((match,val,pslq_res, precision))


        return pslq_results

    def _refine_results(self, intermediate_results, print_results=True):
        return intermediate_results
