"""
Jim DelProposto
Bioinformatics Specialization
optional work.
Code was too slow for a full genome search.
Optimization of previous functions.
Finding Hidden Messages in DNA (Bioinformatics I)
"""

from typing import OrderedDict
import numpy as np
import time
import cProfile
import sys
import os
sys.path.append(os.path.normpath(r"\\wsl.localhost\Ubuntu-20.04\home\delpropo\github\courses\bifx"))
from week1 import extract_seq_from_txt, values_in_range, dl_save_txt


def seq_size_generator(sequence: str, size: int):
    """
    creates a generator which produces all the sequences of the specified size
    This eliminates the need to repeatedly parse through sequences

    param: sequence:  string with the sequence being searched
    param: size:  size of the sequence
    yield: creates a generator which yields all sequences of a specific size
    """
    for i in range(len(sequence) - size + 1):
        yield sequence[i: i + size]

def pattern_locations_np(sequence: str, pattern: str) -> list:
    """
    Make faster using an np array
    searches a sequences for the starting bp locations
    location start with bp as 0, not 1

    param: sequence:  string with the sequence being searched
    param: pattern:  the sequences being searched
    return: a list of the bp locations of the start of the sequence
    >>> " ".join([str(j) for j in pattern_locations("GATATATGCATATACTT", "ATAT")])
    '1 3 9'
    >>> " ".join([str(j) for j in pattern_locations("ATGCATGC", "ATGC")])
    '0 4'
    >>> " ".join([str(j) for j in pattern_locations("ATATATAT", "ATAT")])
    '0 2 4'
    """
    sequence = np.array(list(sequence))
    pattern_len = len(pattern)
    pattern = np.array(list(pattern))
    locations = []

    for i in range(len(sequence) - pattern_len + 1):
        if np.all(pattern == sequence[i:(i + pattern_len)]):
            locations.append(i)
    return locations


def word_starts(sequence: str, size: int) -> dict:
    """
    search for sequences of specific sizes and populate a dictionary with a list of start sites
    returns a dictionary with the locations where the sequence starts
    param: sequence:  string with the sequence being searched
    param: size:  size of the sequences being searched for
    return: a dictionary with a list of the start sequences bases


    """
    # initialize a dictionary to contains the k-mer and sequence start
    seq_starts_dict = {}

    bp = 0
    for i in seq_size_generator(sequence, size):
        # all values are at the same size
        assert len(i) == size
        if i in seq_starts_dict.keys():
            seq_starts_dict[i].append(bp)
        elif i is not None:
            seq_starts_dict[i] =[bp]
        bp += 1

    return seq_starts_dict



def values_in_range_complete(value_list: list, bp_range: int, pattern: str) -> int:
    """
    Determines the maximum number of values in a specified range of a list
    [0, 1, 2, 5, 10] with a range of 2 would have 2 values

    param: value_list: list of values which should be in numerical order
    param: bp_range the largest difference between two values in the list allowed
    return:  max_value:  the largest number of values in a slice where the starting and final value are less than or equal to range_num
    >>> values_in_range_complete([0, 10, 50, 85, 86, 91], 5)
    2
    >>> values_in_range_complete([0, 10, 50, 84, 85, 86, 91], 6, 2)
    3
    >>> values_in_range_complete([0, 1, 2, 3, 4, 5, 6], 6, 10)
    6
    >>> values_in_range_complete([1, 2, 3, 4, 5, 6], 6)
    6
    >>> values_in_range_complete([0, 499, 500], 6)
    2
    """
    #sort(values_in_range)
    num_values = len(value_list)
    length_seq = len(pattern)
    # create a default of 0 for the maximum number of starts is a given bp_range
    max_value = 0
    for i in range(num_values):
        for j in range(i+1, num_values):
            assert j > i and value_list[j] > value_list[i]
            if bp_range >= (value_list[j] - value_list[i] + length_seq) and max_value < j-i + 1:
                max_value = j-i + 1

    return max_value




def main():
    """
    seq_test = "ATGCGCATGCACATGCGCAGTCGTGTGAGATCAGATATGCGCAGCGATGCACGTCGTGTGAATGCGCAGATCAGATGCGATGCACGTCGTGTGAGATCAGATGCGATGCACGTCGTGTGAGATCAGATGCGATGCACGTCGTGTGAGATCAGATGCGATGCACGTCGTGTGAGATCAGATGCGATGCACGTCGTGTGAGATCAGATGCGCA"
    pattern = "ATGCGCA"

    # test clump search on ecoli genome
    url_ecoli = "http://bioinformaticsalgorithms.com/data/realdatasets/Rearrangements/E_coli.txt"
    ecoli_data_path = dl_save_txt(url_ecoli, folder=r".\data\clump_search")
    ecoli_genome = extract_seq_from_txt(ecoli_data_path).split("\n")[0]
    #print(search_clumps(ecoli_genome, 9, 500, 3))

    #print(len(cholera_genome))
    #print((len(ecoli_genome)))
    #print(timeit.timeit('len("cholera_genome")'))
    #print(time.time_ns.split("\n")
    #search_clumps(cholera_genome, 9, 500, 3)
    #print(max_freq_words(cholera_genome, 9, all_values=True))
    #print(search_clumps(cholera_genome, 9, 500, 3))
    print("Starting ECOLI CLUMP SEARCH")
    #print(search_clumps_np(cholera_genome, 9, 500, 3))
    #search_clumps_np(cholera_genome, 9, 500, 3)

    #print(pattern_locations_np(seq_test, pattern))

    cholera_file = os.path.normpath(".\data\PatternMatching_locations\inputs\Vibrio_cholerae.txt")
    print(cholera_file)
    cholera_genome = extract_seq_from_txt(cholera_file).split("\n")[0]
    print(len(ecoli_genome))
    """
    #dict_starts = word_starts(ecoli_genome, 9)
    print("most frequenct 3mere")
    print(word_starts("AAACATAGGATCAAC", 2))

    #with open(file_path, "w") as text_file:
    #    text_file.write(sequence)

    """
    count_seq = 0
    total_seq = 0
    for key, value in dict_starts.items():
        num_in_range = values_in_range(value, 500)
        total_seq += 1
        if num_in_range > 2:
            count_seq += 1




    print(count_seq)
    print(total_seq)
    #with open(file_path, "w") as text_file:
    #    text_file.write(sequence)

    print("correct range")
    correct_range_count = 0
    for key, value in dict_starts.items():
        if values_in_range_complete(value, 500, key) > 2:
            correct_range_count += 1

    print(correct_range_count)

"""
if __name__ == "__main__":
    main()
    #import doctest
    #doctest.testmod

