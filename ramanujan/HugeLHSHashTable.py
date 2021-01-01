from .LHSHashTable import *
import csv
import re
import os
import pandas as pd

class HugeLHSHashTable(LHSHashTable):
    """
    This class modifies the way LHSHashTable saves all lhs_possibilities, to allow it to be huge
    (bigger then what is possible to store in RAM)

    it is done by splitting the domain to partitions, each partition is stored in different files,
    and has a specific range of values that exp_value might match to (*)
    partition might be splitted to parts, without any specific order between each part.
    This means we can load only one partition when looking for values, and hold less data in the 
    memory.

    This logic is complemented by a matching enumerator that queries this class with sorted items,
    allowing us to read every partition only once, reducing disk calls.


    for example
        zeta.0.0.csv
        zeta.0.1.csv
        zeta.1.0.csv
        zeta.1.1.csv
    there are two partitions (the first number), and each partition has 2 parts.

    (*) the there is no min or max value that might be accepted, so we will sort according to 
        the first 10 digits after the dot.
    """
    def __init__(self, name, search_range, const_vals, key_digits, 
        partitions=30, max_mem_useage_mb=512, force_reload=False, 
        file_name_format='{folder}/{partition}.{chunk}.csv',
        skip_bloom_load=False) -> None:
        """
        there is no good reason not to change father init and allow easy
        overriding, just lazy

        new params:
            partitions - number of partitions to split to. the more partitions there are
                more files and disk calls, but searching each partition will be faster
            max_mem_useage_mb - according to this, decide on the maximum part size
            force_reload - ignore stored values, and enumerate the domain again
        """
        self.name = name
        self.file_name_format = file_name_format
        self.threshold = 10 ** (-key_digits)
        self.key_factor = 10 ** key_digits
        self.max_key_length = key_digits * 2
        constants = [mpmath.mpf(1)] + const_vals
        self.n_constants = len(constants)
        self.partitions = partitions
        self.search_range = search_range

        self.max_capacity = (search_range * 2 + 1) ** (self.n_constants * 2)
        self.bloom = BloomFilter(capacity=self.max_capacity, error_rate=0.05)
        
        start_time = time()

        # calculate the approximate size of each entry, and split to partition. 
        # each has one chunk loaded in mem while enumerating
        self.max_chunk_entries = int((max_mem_useage_mb * 1024 * 1024 / 30) / partitions)

        # the number of unique items possible in each partition
        self.partition_max_unique_ids = int(self.key_factor / partitions)

        self.loaded_partition = pd.DataFrame(columns=['coef_top', 'coef_bottom'])
        self.loaded_partition_id = -1

        self.partition_chunk_map = {}
        # not a good test
        if force_reload or not os.path.isfile(file_name_format.format(
            folder=self.name, partition=0, chunk=0)):
            self._init_lhs_items_cache()
            self._enumerate_lhs_domain(constants, search_range, self.key_factor)
        else:
            self._load_from_folder(skip_bloom_load)
            
        print('initializing LHS dict: {}'.format(time() - start_time))

    def _load_from_folder(self, skip_bloom_load):
        # we want to extract the partition id from the name
        name_format_no_folder = self.file_name_format.split('/')[1]
        re_part = name_format_no_folder.replace('{chunk}','\d+').replace('{partition}','(\d+)')

        for file in os.listdir(self.name):
            try:
                partition_id = re.findall(re_part, file)[0]
                self._add_to_partition_chunk_map(file, partition_id)

                if skip_bloom_load:
                    continue

                with open(self.name + '/' + file, 'r') as f:
                    csv_reader = csv.reader(f)
                    for line in csv_reader:
                        self.bloom.add(line[0])


            except Exception as e:
                print(f'Error on loading from {self.name}/{file}.\n{e}\ncontinuing...')

    def _enumerate_lhs_domain(self, *args, **kwargs):
        print('creating domain')
        super()._enumerate_lhs_domain(*args, **kwargs)

        # since super()._enumerate_lhs_domain assumes that self._add_to_lhs_possblilites
        # commits data instantly, and under this subclass, the data is committed when a 
        # threshold is passed, the last chunk might never be committed. Committing explicitly
        for i in range(self.partitions):
            if len(self.lhs_items_cache[i]) != 0:
                self._flush_partition_file(i)

        # at this point, self.lhs_items_cache should be empty. This will reduce memory usage

    def _get_by_key(self, key):
        # create a o(1) formula that finds the file id, and open it and only it
        # consider not dumping the data, but storing it, assuming that the next calls are likely
        # to be in the same file
        partition_id = self._get_partition_id(str(key))
        if partition_id != self.loaded_partition_id:
            # get next partition
            print(f'loading partition no. {partition_id}')
            self.loaded_partition = pd.DataFrame(columns=['coef_top', 'coef_bottom'])
            self.loaded_partition_id = partition_id

            for chunk_file_name in self.partition_chunk_map[partition_id]:
                self._load_csv_to_cache(self.name + '/' + chunk_file_name)
                
        lines = self.loaded_partition.loc[int(key)].values
        matches = []
        try:
            # fix me
            c_top_str, c_bottom_str = lines
            c_top = tuple([int(i) for i in c_top_str.split(',')])
            c_bottom = tuple([int(i) for i in c_bottom_str.split(',')])
            return [(c_top, c_bottom)]
        
        except Exception as e:
            print(e)
            print(key)
            print(lines)

            for c_top_str, c_bottom_str in self.loaded_partition.loc[int(key)].values:
                c_top = tuple([int(i) for i in c_top_str.split(',')])
                c_bottom = tuple([int(i) for i in c_bottom_str.split(',')])

                matches.append((c_top, c_bottom))

        return matches

    def _load_csv_to_cache(self, csv_path):
        chunk = pd.read_csv(csv_path, names=['coef_top', 'coef_bottom'])
        self.loaded_partition = self.loaded_partition.append(chunk)
        # with open(csv_path, 'r') as f:
        #     csv_reader = csv.reader(f)
        #     for exp_value, c_top_str, c_bottom_str in csv_reader:
        #         exp_value = int(exp_value)
        #         c_top = tuple([int(i) for i in c_top_str.split(',')])
        #         c_bottom = tuple([int(i) for i in c_bottom_str.split(',')])

        #         compressed_coefs = self._compres_coefs(c_top, c_bottom)

        #         if exp_value in self.loaded_partition:
        #             self.loaded_partition[exp_value].append(compressed_coefs)
        #         else:
        #             self.loaded_partition[exp_value] = [compressed_coefs]

    def _add_to_lhs_possblilites(self, exp_value, c_top, c_bottom):
        # get name and subsitude to name_format
        partition_id = self._get_partition_id(exp_value)
        c_top_str = ','.join([str(i) for i in c_top])
        c_bottom_str = ','.join([str(i) for i in c_bottom])

        self.lhs_items_cache[partition_id].append((exp_value, c_top_str, c_bottom_str))

        if len(self.lhs_items_cache[partition_id]) == self.max_chunk_entries:
            self._flush_partition_file(partition_id)

    def _flush_partition_file(self, partition_id):
        # keep different files that hold the data by key ranges
        print(f'flushing {partition_id}. \n {[len(self.lhs_items_cache[i]) for i in self.lhs_items_cache]}')
        chunk_file_name = self.file_name_format.format(folder=self.name, partition=partition_id, 
            chunk=self.partition_currrent_chunk[partition_id])

        with open(chunk_file_name, 'w') as f:
            writer = csv.writer(f)
            writer.writerows(self.lhs_items_cache[partition_id])

            # layout for next insertion
            self.lhs_items_cache[partition_id] = []
            self.partition_currrent_chunk[partition_id] += 1

            self._add_to_partition_chunk_map(chunk_file_name, partition_id)

    def _get_partition_id(self, exp_str):
        exp_str = exp_str if exp_str[0] != '-' else exp_str[1:]
        exp_sort_val = int(exp_str[-len(str(int(self.key_factor)))+1:])

        return int(exp_sort_val / self.partition_max_unique_ids)

    def _add_to_partition_chunk_map(self, chunk_file_name, partition_id):
        partition_id = int(partition_id)
        if partition_id in self.partition_chunk_map:
            self.partition_chunk_map[partition_id].append(chunk_file_name)
        else:
            self.partition_chunk_map[partition_id] = [chunk_file_name]

    def _init_lhs_items_cache(self):
        try:
            os.mkdir(self.name)
        except FileExistsError:
            pass 

        self.lhs_items_cache = {}
        for i in range(self.partitions):
            self.lhs_items_cache[i] = []

        self.partition_currrent_chunk = [0 for i in range(self.partitions)]