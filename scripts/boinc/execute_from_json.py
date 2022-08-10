import os
import sys
import json
from ramanujan.poly_domains.Zeta3Domain1 import Zeta3Domain1
from ramanujan.poly_domains.Zeta3Domain2 import Zeta3Domain2
from ramanujan.poly_domains.Zeta5Domain import Zeta5Domain
from ramanujan.poly_domains.Zeta7Domain import Zeta7Domain
from ramanujan.poly_domains.CatalanDomain import CatalanDomain
from ramanujan.poly_domains.ExplicitCartesianProductPolyDomain import (
    ExplicitCartesianProductPolyDomain,
)
from ramanujan.constants import g_const_dict
from ramanujan.enumerators.FREnumerator import FREnumerator


ENUMERATORS = {"FREnumerator": FREnumerator}
DOMAINS = {
    "Zeta3Domain1": Zeta3Domain1,
    "Zeta3Domain2": Zeta3Domain2,
    "Zeta5Domain": Zeta5Domain,
    "Zeta7Domain": Zeta7Domain,
    "CatalanDomain": CatalanDomain,
    "ExplicitCartesianProductPolyDomain": ExplicitCartesianProductPolyDomain,
}
RESULTS_FILE_PATH_FORMAT = "{id}_results.json"


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
    unique_id = os.path.split(json_path)[1].rsplit(".", 1)[0]

    with open(json_path, "r") as f:
        config = json.load(f)

    # some domains accept their ranges in different ways.
    # (If the scheme only contains one coef range, it will be on list. If it has more ranges, it may be a list of lists)
    # This way we pass the non-standard init function and assign the ranges directly
    poly_domain = DOMAINS[config["domain_type"]]()
    poly_domain.a_coef_range = config["an_coefs"]
    poly_domain.b_coef_range = config["bn_coefs"]
    poly_domain.only_balanced_degrees = config["only_balanced_degrees"]
    poly_domain.use_strict_convergence_cond = config["use_strict_convergence_cond"]
    poly_domain._setup_metadata()

    const_vals = get_consts_objects(config["const_list"])
    enumerator = ENUMERATORS[config["enumerator"]](poly_domain, const_vals)

    results = enumerator.find_initial_hits()

    with open(RESULTS_FILE_PATH_FORMAT.format(id=unique_id), "w") as f:
        json.dump(results, f)


if __name__ == "__main__":
    main()
