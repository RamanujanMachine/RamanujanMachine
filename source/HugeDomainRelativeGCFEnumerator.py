from RelativeGCFEnumerator import *
from HugeLHSHashTable import *

class HugeDomainRelativeGCFEnumerator(RelativeGCFEnumerator):
    """
        when fiddling with huge domains, enumerating the RHS options, using 
        LHSHashTable might exceed process max memory. This class uses the 
        HugeLHSHashTable, that handles the heavy lifting, and making minor
        changes to RelativeGCFEnumerator to match it.
    """

    def __init__(self, *args, **kwargs):
        print('using relative enumerator')
        super().__init__(*args, **kwargs)
        
    def _init_hash_table(self, saved_hash):
        with mpmath.workdps(self.enum_dps):
            constants = [const() for const in self.constants_generator]
            self.hash_table = HugeLHSHashTable(
                saved_hash, 
                self.lhs_limit,
                constants, 
                g_N_initial_key_length)

    def _refine_results(self, intermediate_results: List[Match], print_results=True):
        intermediate_results.sort(key=lambda x: self.hash_table._get_partition_id(str(x.lhs_key)))
        return super()._refine_results(intermediate_results, print_results)