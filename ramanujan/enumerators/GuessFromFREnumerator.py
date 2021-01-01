import math
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl, wlexpr

from .FREnumerator import *


class GuessFromFREnumerator(FREnumerator):
    """
        Uses FREnumerator to get all GCF that have FR, 
        Then uses wolfarm-api to try and evaluate the expression 
        and get a closed form from it
    """

    def __init__(self, *args, **kwargs):
        print('Using wolfarm to eval expressions')
        self.session = WolframLanguageSession(
            '/home/jellybean/wolfarm/Executables/WolframKernel')
        super().__init__(*args, **kwargs)
    

    def _refine_results(self, intermediate_results: List[Match], print_results=True):
        results = []

        for res in intermediate_results:
            # calculate the gcf a bit better
            an_iter_func, bn_iter_func = self.poly_domains_generator.get_calculation_method()
            an_iter = an_iter_func(res.an_coef, 3000, start_n=0)
            bn_iter = bn_iter_func(res.bn_coef, 3000, start_n=0)

            gcf = gcf_calculation_to_precision(an_iter, bn_iter, 20, 
                min_iters=2000, burst_number=500)
            rhs_str = mpmath.nstr(gcf, 20)
            
            results.append(
                self.session.evaluate(
                    wlexpr(
                        f'WolframAlpha["{rhs_str}", IncludePods -> "PossibleClosedForm"]'
                        )
                    )
                )

        return results