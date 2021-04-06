from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.constants import g_const_dict
from ramanujan.multiprocess_enumeration import multiprocess_enumeration

"""
This script shows an usage example for multiprocessing.
It mimics the file paper_results/zeta3_results, but using multiprocessing.

Please address the docs under paper_results/zeta3_results for more information.
"""


def main():
    # creating the LHS is done the same way.
    saved_hash = 'paper_results/zeta3.lhs.dept20.db'
    lhs_search_limit = 20
    lhs = LHSHashTable(
        saved_hash,
        lhs_search_limit,
        [g_const_dict['zeta'](3)])

    # You can now define a much larger domain, and get results in a decent time.
    # The domain definition stays the same. We will split it later to chunks and enumerate each chink
    # in a different process
    poly_search_domain = Zeta3Domain1(
        [(2, 2), (1, 1), (-100, 100), (-100, 100)],
        (-50, -1))

    # Instead of creating the enumerator yourself, pass the type of enumerator you wish to use, and
    # all of the other parameters you need to create an execution. Specify the number of processes to 
    # spawn. We will create split the domains to chunks according to the required number of processes,
    # and create an enumerator for each chunk. 
    # the following example will spawn 4 process, using EfficientGCFEnumerator
    results = multiprocess_enumeration(
        EfficientGCFEnumerator,
        lhs,
        poly_search_domain,
        [g_const_dict['pi']],
        4)

    print(results)


if __name__ == '__main__':
    main()
