import multiprocessing
from copy import deepcopy


class Dummy(object):
    """
    When passing the lhs object to each process, only the bloom filter is used.
    To pass less data to each process, we'll be using this class to duplicate only
    the bloom filter. 
    Also, passing LHSHashTable to a child process requires us to pickle it, which creates
    problems. This way, we only handle the required data when spawnning a new process
    """
    pass


def _single_process_execution(enumerator_class, lhs, poly_search_domain, const_vals):
    import os
    print(os.getpid())
    enumerator = enumerator_class(
        lhs,
        poly_search_domain,
        const_vals)
    
    return enumerator.find_initial_hits()

def multiprocess_enumeration(enumerator_class, lhs, poly_search_domain, const_vals, number_of_processes):
    """
    This function will split an execution to number_of_processes different processes
    The poly_domain will be splitted, and for each chunk an instance of lhs and 
    enumerator will be created. Each instance will preform the first enumeration.
    The refinning process requires to load the LHS dict to memory. This will be done by
    only one instance, when all first enumerations are finnished.

    :param enumerator_class: the CLASS (NOT an instance) of the requested enumerator
    :param lhs: an LHSHashTable object
    :param poly_search_domain: the requested PolyDomain - this object will be splitted to processes
    """

    pool = multiprocessing.Pool(processes=number_of_processes)
    arguments = []
    
    # Each subprocess only uses lhs.bloom. See Dummy class doc for more details.
    lean_lhs = Dummy()
    lean_lhs.bloom = lhs.bloom

    # Creating arguments for each process function
    splitted_domain = poly_search_domain.split_domains_to_processes(number_of_processes)
    for domain_chunk in splitted_domain:
        arguments.append((
            enumerator_class,
            deepcopy(lean_lhs.bloom),
            domain_chunk,
            const_vals
            ))

    process_results = pool.starmap(_single_process_execution, arguments)
    unified_results = []
    for r in process_results:
        unified_results += r

    # Create another enumerator (should not take time to initiate) and preforme 
    # the second step using only it
    enumerator = enumerator_class(lhs, poly_search_domain, const_vals)
    refined_results = enumerator.refine_results(unified_results)

    return refined_results