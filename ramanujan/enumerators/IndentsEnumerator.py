from .FREnumerator import FREnumerator, Match, RefinedMatch
import mpmath
from collections import namedtuple
from .RelativeGCFEnumerator import RelativeGCFEnumerator


CONST_LIST = ['zeta(3)', 'zeta(4)', 'zeta(5)', 'pi', 'pi*pi', 'e']

RefinedMatch2 = namedtuple('RefinedMatch2', 'rhs_an_poly rhs_bn_poly val expr precision')


class IndentsEnumerator(FREnumerator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _first_enumeration(self, print_results: bool):
        results = []  # list of intermediate results        
        for an_iter, bn_iter, metadata in self._iter_domains_with_cache(1):
            # everyone have fr
            results.append(Match(metadata.an_coef, metadata.bn_coef))
        return results

    def _identify_results(self, precise_intermediate_results, const_list=CONST_LIST, verbose=True):
        results = []
        for match, val, precision in precise_intermediate_results:
            mpf_val = mpmath.mpf(val)
            print(match)
            try:
                identified_res = mpmath.identify(mpf_val, const_list)
            except Exception as e:
                import ipdb
                ipdb.set_trace()
            if identified_res:
                print(identified_res)
                results.append(RefinedMatch(*match, val, identified_res, precision))
            else:
                results.append(RefinedMatch(*match, val, None, precision))

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
        pslq_res = super()._improve_results_precision(intermediate_results, verbose)

        print('Trying to identify results')
        identify_res = []
        for match, val, precision in self.precise_intermediate_results:
            mpf_val = mpmath.mpf(val)
            print(match)
            try:
                identified_res = mpmath.identify(mpf_val, const_list)
            except Exception as e:
                identified_res = False
                # import ipdb
                # ipdb.set_trace()
            if identified_res:
                print(identified_res)
                identify_res.append(RefinedMatch2(*match, val, identified_res, precision))
            else:
                identify_res.append(RefinedMatch2(*match, val, None, precision))

        return pslq_res + identify_res

    def _refine_results(self, intermediate_results, print_results=True):
        return intermediate_results
