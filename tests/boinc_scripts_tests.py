import unittest
import json
import sys
import os
from ramanujan.poly_domains.Zeta5Domain import Zeta5Domain
from unittest.mock import patch
from shutil import rmtree

# Since the BOINC scripts are not a part of the ramanujan package, we'll add their folder to python's path.
# We must do this step before we can import the functions we wish to test, since we do not export an API to change
# those script's configurations - its defined in their main function.
# This can cause IDEs to declare some lines as an unresolved reference. Also, this is  breaking coding conventions :(
current_path = os.getcwd()
boinc_scripts_dir = os.path.join(current_path, '..', 'scripts', 'boinc')
sys.path.insert(1, boinc_scripts_dir)
from split_execution import split_to_jsons, store_execution_to_json
import execute_from_json


class BoincTests(unittest.TestCase):
    def test_split_execution(self):
        """
        Test the Split to files option in BOINC.
        This includes comparing the contents of each file, and making sure they contains the entire original domain,
        and that parameters are correct
        """
        const_list = [['zeta', 5], ['zeta', 3]]
        identifier = 'boinc_script_tests'
        # delete folder from previous tests
        if os.path.isdir(identifier):
            rmtree(identifier)

        poly_search_domain = Zeta5Domain(
            [(1, 1), (-100, 100), (-100, 100)],
            (1, 1),
            only_balanced_degrees=True)

        split_to_jsons(identifier, "FREnumerator", poly_search_domain, const_list)

        split_ranges = []
        for config_file in os.listdir(identifier):
            with open(os.path.join(identifier, config_file), 'r') as f:
                config = json.load(f)
            self.assertEqual(config['enumerator'], "FREnumerator")
            self.assertEqual(config['only_balanced_degrees'], True)
            self.assertEqual(config['use_strict_convergence_cond'], False)
            self.assertEqual(config['const_list'], const_list)
            split_ranges.append([config['an_coefs'], config['bn_coefs']])

        # make sure all items are included in the split domains
        original_items = [i for i in poly_search_domain.iter_polys('a')]
        split_items = []
        for i in split_ranges:
            split_domain = Zeta5Domain(
                i[0], i[1][0],
                only_balanced_degrees=True)
            split_items += [i for i in split_domain.iter_polys('a')]
        original_items.sort()
        split_items.sort()
        self.assertEqual(split_items, original_items)

    def test_load_execute_from_json(self):
        example_json_filename = 'boinc_example_config.json'
        expected_result_filename = 'boinc_example_config_results.json'
        # remove files from previous executions
        if os.path.isfile(expected_result_filename):
            os.remove(expected_result_filename)

        with patch.object(sys, 'argv', ['execute_from_json.py', example_json_filename]):
            execute_from_json.main()

        self.assertTrue(os.path.isfile('boinc_example_config_results.json'))
        with open(expected_result_filename, 'r') as f:
            results = json.load(f)
        self.assertIn(
            [[5, 4, 1], [-4, 2, 0, 0, 0],
             "0.7478010885109415849488278024886662133297666830275923650467536636688412927137336490450574799505194538",
             [], [], 100],
            results)
        self.assertIn(
            [[5, 6, 2], [-4, 2, 0, 0, 0],
             "1.823781305562079885989830337775097500278457944100437879220962574095116177320237941441766930661193077",
             [18, 0], [0, 1], 100],
            results)


if __name__ == '__main__':
    unittest.main()
