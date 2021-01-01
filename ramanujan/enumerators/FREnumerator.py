import math
from .RelativeGCFEnumerator import *

CONVERGENCE_THRESHOLD = 0.1

def check_for_fr(an_iterator, bn_iterator, metadata, burst_number=200, min_iters=1):
    deg_a = len(metadata.an_coef) - 1
    deg_b = len(metadata.bn_coef) - 1
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
                    (deg_a) * (-mpmath.log(i) + 1)
            )

            if num_of_calculated_vals >= 3 and \
                abs(calculated_values[-2] - calculated_values[-1]) > abs(calculated_values[-2] - calculated_values[-3]):
                return False, i

            if num_of_calculated_vals >= 2 and \
                abs(calculated_values[-2] - calculated_values[-1]) < CONVERGENCE_THRESHOLD:
                return True , i


    return False, i


class FREnumerator(RelativeGCFEnumerator):
    """
        This enumerator calculates the gcf to an arbitrary dept, until getting to
        a stable value withing the precession required bounds.
        Useful for GCF that converges slowly.

        Note-
            Currently using costume series generators is not allowed, since those 
            operate by computing a fixed n items, as opposed to the fluid approach
            which suggested here. To implement it, one need to make those an iterative
            function, and pass it to cached_series_iterator (and also add support there)
    """

    def __init__(self, *args, **kwargs):
        print('checking for FR enumerator')
        super().__init__(*args, **kwargs)


    def _first_enumeration(self, poly_domains, print_results: bool):
        num_iterations = poly_domains.num_iterations

        start = time()
        key_factor = 1 / self.threshold

        results = []  # list of intermediate results        
        all_items_calculated = []
        last_lead = -1001 # random number that will not pop
        for an_iter, bn_iter, metadata in self._gcf_series_cached_iters(poly_domains, 100000):
            if metadata.an_coef[0] != last_lead:
                last_lead = metadata.an_coef[0]
                print(f"last_lead {last_lead} - {time() - start}")
            has_fr, items_calculated = check_for_fr(an_iter, bn_iter, metadata)
            if has_fr:
                all_items_calculated.append(items_calculated)
                print (metadata)
                results.append(metadata)

        if len(all_items_calculated) != 0:
            print(sum(all_items_calculated)/len(all_items_calculated))


        return results

    def _refine_results(self, intermediate_results: List[Match], print_results=True):
        pass