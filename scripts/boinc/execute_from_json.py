import json
import os
import sys
import json
import time
import ramanujan
import pybloom_live
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.constants import g_const_dict
from ramanujan.LHSHashTable import LHSHashTable
from ramanujan.enumerators.EfficientGCFEnumerator import EfficientGCFEnumerator

# we're gonna have to load all of the domains :(
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1

enumerators = {
    "EfficientGCFEnumerator": EfficientGCFEnumerator
}
domains = {
    "Zeta3Domain1": Zeta3Domain1
}

class Dummy(object):
    pass

def get_consts_objects(consts):
    const_objects = []
    for i in consts:
        if isinstance(i, tuple) or isinstance(i, list):
            const_objects.append(g_const_dict[i[0]](i[1]))
        else:
            const_objects.append(g_const_dict[i])
    return const_objects

def main():
    json_path = sys.argv[1]
    bloom_path = sys.argv[2]

    with open(json_path, 'r') as f:
        config = json.load(f)

    lean_lhs = Dummy()
    with open(bloom_path, 'rb') as f:
        lean_lhs.bloom = pybloom_live.BloomFilter.fromfile(f)

    poly_domain = domains[config["domain_type"]]()
    poly_domain.a_coef_range = config["an_coefs"]
    poly_domain.b_coef_range = config["bn_coefs"]
    poly_domain._setup_metadata()

    const_vals = get_consts_objects(config["const_list"])

    enumerator = enumerators[config["enumerator"]](
        lean_lhs.bloom,
        poly_domain,
        const_vals)

    results = enumerator.find_initial_hits()
        print(results)

    # TODO - write results to file

if __name__ == '__main__':
    main()