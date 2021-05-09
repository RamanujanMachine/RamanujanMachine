import os
import json
import time
import ramanujan
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.constants import g_const_dict
from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1


"""
In this script, we create multiple jobs from a single execution scheme. 
For each job save a json that stores the execution parameters. 
Those jsons will be sent to different hosts from a BOINC server

The bloom filter is saved 
"""
SPLITTED_DOMAIN_CHUNK_SIZE = 100_000
BLOOM_LOCAL_PATH = "./{0}_bloom.bin"
JSON_FORMAT = "./{folder}/{id}.json"

ALLOWD_ENUMERATORS = [
    "EfficientGCFEnumerator"
]

def dump_bloom_to_file(identifier, lhs):
    with open(BLOOM_LOCAL_PATH.format(identifier), 'wb') as f:
        lhs.bloom.tofile(f)


def to_json(identifier, enumerator_type, domain, domain_type, const_list):
    number_of_chunks = poly_search_domain.num_iterations // SPLITTED_DOMAIN_CHUNK_SIZE
    domain_chunks = poly_search_domain.split_domains_to_processes(number_of_chunks)
    id_max_len = len(str(len(domain_chunks)))
    os.mkdir(identifier)

    for i, chunk in enumerate(domain_chunks):
        dest_file_name = JSON_FORMAT.format(
            folder=identifier, 
            id=str(i).zfill(id_max_len))
        with open(dest_file_name, 'w') as f: 
            chunk_data = {
                "an_coefs": chunk.a_coef_range,
                "bn_coefs": chunk.b_coef_range,
                "enumerator": enumerator_type,
                "domain_type": domain_type,
                "const_list": const_list
            }
            json.dump(chunk_data,f)
            
def get_consts_objects(consts):
    const_objects = []
    for i in consts:
        if isinstance(i, tuple) or isinstance(i, list):
            const_objects.append(g_const_dict[i[0]](i[1]))
        else:
            const_objects.append(g_const_dict[i])
    return const_objects

    
def main():
    # Change the following arguments as you wish 
    const_list = ['e']
    saved_hash = 'e_lhs_dept5_db'
    lhs_search_limit = 5
    identifier = None # distinct name for generated json

    const_objects = get_consts_objects(const_list)
    identifier = identifier if identifier else time.strftime("%y%m%d_%H%M%S")
    lhs = LHSHashTable(
        saved_hash,
        lhs_search_limit,
        const_objects) 

    poly_search_domain = Zeta3Domain1(
        [(2, 2), (1, 1), (1, 100), (1, 100)],
        (-50, -1))

    to_json(identifier, "EfficientGCFEnumerator", poly_search_domain, "Zeta3Domain1", const_list)
    EfficientGCFEnumerator,
        lhs,
        poly_search_domain,
        const_objects


if __name__ == '__main__':
    main()