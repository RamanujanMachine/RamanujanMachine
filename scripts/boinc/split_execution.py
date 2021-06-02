import os
import json
import time
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.poly_domains.Zeta3Domain2 import Zeta3Domain2


"""
In this script, we create multiple jobs from a single execution scheme. 
For each job save a json that stores the execution parameters. 
Those jsons will be sent to different hosts from a BOINC server

Currently, we're planning on distributing FREnumerator only, so there is no support for 
bloom filter distribution.
"""
SPLIT_DOMAIN_CHUNK_SIZE = 100
BLOOM_LOCAL_PATH = "./{0}_bloom.bin"
JSON_NAME_FORMAT = "./{folder}/{filename}.json"

ALLOWED_ENUMERATORS = [
    "FREnumerator"
]
POLY_DOMAINS_MAPPING = {
    'Zeta3Domain1': Zeta3Domain1,
    'Zeta3Domain2': Zeta3Domain2
}


def store_execution_to_json(dest_file_name, enumerator_type, poly_domain, const_list):
    with open(dest_file_name, 'w') as f: 
        chunk_data = {
            "an_coefs": poly_domain.a_coef_range,
            "bn_coefs": poly_domain.b_coef_range,
            "enumerator": enumerator_type,
            "domain_type": poly_domain.__class__.__name__,
            "const_list": const_list
        }
        json.dump(chunk_data, f)


def split_to_jsons(identifier, enumerator_type, domain, const_list):
    if enumerator_type not in ALLOWED_ENUMERATORS:
        raise ValueError(f"Required enumerator is not supported. Use one from the following:\n{ALLOWED_ENUMERATORS}")
    number_of_chunks = domain.num_iterations // SPLIT_DOMAIN_CHUNK_SIZE
    domain_chunks = domain.split_domains_to_processes(number_of_chunks)
    id_max_len = len(str(len(domain_chunks)))
    os.mkdir(identifier)

    for i, chunk in enumerate(domain_chunks):
        dest_file_name = JSON_NAME_FORMAT.format(
            folder=identifier, 
            filename=f"{identifier}_{str(i).zfill(id_max_len)}")
        store_execution_to_json(dest_file_name, enumerator_type, chunk, const_list)


def main():
    # Change the following arguments as you wish 
    const_list = [('zeta', 3)]
    identifier = None  # distinct name for generated json

    identifier = identifier if identifier else time.strftime("%y%m%d_%H%M%S")
    
    poly_search_domain = Zeta3Domain2(
        [(1, 100), (1, 100)],
        (1, 50))

    split_to_jsons(identifier, "FREnumerator", poly_search_domain, const_list)


if __name__ == '__main__':
    main()
