

import hydra
import sys
import os

sys.path.append(os.path.normpath(r"\\wsl.localhost\Ubuntu-20.04\home\delpropo\github\courses\bifx"))
from week1 import extract_seq_from_txt, max_values_in_range, dl_save_txt, DNA_BASES
#from genome_search_optional import

SKEW_DICT = {
            "C": -1,
            "G":  1,
            "A":  0,
            "T":  0,
            }


def folder_txt_files(folder: str) -> tuple:
    """
    returns a list of paths of all text files in the specified folder
    :param str name: folder path
    :returns: list of files in the folder that end with .txt

    """
    path = os.path.normpath(folder)
    file_list = [os.path.join(path, file) for file in os.listdir(path) if file.endswith(".txt")]
    return (file_list)



def answer_print(input: str) -> str:
    """
    joins values in a s

    """
    print("".join(str(x) + " " for x in input).strip())



def seq_GC_skew(sequence: str) -> str:
    """
    input a sequence as a string
    return the skew valve for each base in the string
    :param: sequence: DNA sequence
    :return: integers which give the skew for each G or C base

    >>> answer_print(seq_GC_skew("CATGGGCATCGGCCATACGCC"))
    0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
    """



    skew_count = 0
    skew_output = [0]
    for base in sequence:
        if base in SKEW_DICT.keys():
            skew_count += SKEW_DICT[base]
            skew_output.append(skew_count)
    return skew_output


def min_skew_index(skew_values: list) -> list:
    """
    input the seq_GC_skew and return a list with the index of the minimum values

    :param: skew_values: list of integers
    :return: list of the locations of the minimum values


    """

    min_value = min(skew_values)
    return [i for i, x in enumerate(skew_values) if x == min_value]

def hamming_distance(seq_one: str, seq_two: str) -> int:
    """
    return the number of mismatches between two sequences
    This is defined as the Hamming distance
    :param: seq_one: The first sequence
    :param: seq_two: The second sequence
    :return: return the hamming distance integer

    >>> hamming_distance("GGGCCGTTGGT","GGACCGTTGAC" )
    3
    """
    assert len(seq_one) == len(seq_two)

    hamming_distance = 0
    for base, i in enumerate(seq_one):
        if i != seq_two[base]:
            hamming_distance += 1
    return hamming_distance

def aproximate_pattern_match(pattern: str, sequence: str, mismatches: int) -> tuple:
    """
    All starting positions where the pattern appears as a substring of text with at most d mismatches

    #>>>  answer_print(aproximate_pattern_match("ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3))
    #6 7 26 27
    """


    assert mismatches > 0

    locations = []
    pattern_len = len(pattern)
    for i in range(len(sequence)-pattern_len + 1):
        if mismatches >= hamming_distance(pattern, sequence[i:(i + pattern_len)]):
            locations.append(i)
    return locations

def approximate_pattern_count(pattern: str, sequence: str, mismatches: int) -> int:
    """
    count the number of values from the approximate_pattern_match function

    >>> approximate_pattern_count("GAGG", "TTTAGAGCCTTCAGAGG", 2)
    4
    """
    return len(aproximate_pattern_match(pattern, sequence, mismatches))



def immediate_neighbors(pattern: str) -> set:
    """
    input a sequence
    returns all sequences where the hamming distance is 1
    find all related sequence where the hamming distance does not exceed the input value
    iterate through each individual base to create a mismatch 1, 2, 3, etc
    recursively call immediate_neighbors to create the other values
    """

    sequence_set = set()
    for i in range(len(pattern)):
        for base in DNA_BASES:
            if pattern[i] != base:
                temp_seq = list(pattern)
                temp_seq[i] = base
                sequence_set.add("".join(temp_seq))

    for value in sequence_set:
        assert hamming_distance(pattern, value) == 1

    return sequence_set

def neighbors(pattern: str, d: int):
    """
    recursively call immediate neighbors on all iterations of neighbors until
    the hamming distance is less than or equal to d

    """
    pass


def frequent_words_with_mismatches(sequence: str, length: int, mismatches: int) -> list:
    """
    count the number of times a string has an approximate match in the text
    for a given kmer substring, increate the count of every kmer
    that has Hamming distance at most d from pattern
    this is called
    d-neighborhood of pattern

    """

    word_list = []
    for i in range(len(sequence - length + 1)):
        pattern = sequence[i:i + length]



    return word_list


@hydra.main(config_path=".\config", config_name="config")
def main(cfg):

    print(cfg.skew_sample.input)
    print(cfg.skew_sample.output)
    text = seq_GC_skew(cfg.skew_sample.input)
    print("answer")
    print(answer_print(text))
    print("test_problem")
    print(answer_print(seq_GC_skew(cfg.skew_test.input)))

    print(seq_GC_skew(cfg.minimum_skew_sample.input))
    print(min_skew_index(seq_GC_skew(cfg.minimum_skew_sample.input)))

    print(seq_GC_skew(cfg.minimum_skew_sample.input))
    text = min_skew_index(seq_GC_skew(cfg.minimum_skew_sample.input))
    answer_print(text)



    # skew test
    skew_test_seq = extract_seq_from_txt(os.path.normpath(cfg.files.min_skew_test))
    print(min_skew_index(seq_GC_skew(skew_test_seq)))

    # Hamming distance sample test
    hamming_sample_txt = extract_seq_from_txt(os.path.normpath(cfg.files.hamming_distance_sample))
    seq_one, seq_two = hamming_sample_txt.split('\n')
    print(hamming_distance(seq_one, seq_two))

    # hamming distance test
    hamming_sample_txt = extract_seq_from_txt(os.path.normpath(cfg.files.hamming_distance_test))
    seq_one, seq_two = hamming_sample_txt.split('\n')
    print(hamming_distance(seq_one, seq_two))

    # approximate pattern matching problem
    apms_sample_txt = extract_seq_from_txt(os.path.normpath(cfg.files.approximate_pattern_match_sample))
    seq_one, seq_two, mismatch = apms_sample_txt.split('\n')
    print(answer_print(aproximate_pattern_match(seq_one, seq_two, int(mismatch))))


    apms_sample_txt = extract_seq_from_txt(os.path.normpath(cfg.files.approximate_pattern_match_test))
    seq_one, seq_two, mismatch = apms_sample_txt.split('\n')
    print(answer_print(aproximate_pattern_match(seq_one, seq_two, int(mismatch))))

    # frequent approximate pattern matching problem
    seq_one = cfg.freq_approximate_pattern_match_sample.pattern
    seq_two = cfg.freq_approximate_pattern_match_sample.sequence
    mismatches = cfg.freq_approximate_pattern_match_sample.mismatches
    print(len(aproximate_pattern_match(seq_one, seq_two, int(mismatches))))

    seq_one = cfg.freq_approximate_pattern_match_test.pattern
    seq_two = cfg.freq_approximate_pattern_match_test.sequence
    mismatches = cfg.freq_approximate_pattern_match_test.mismatches
    print(len(aproximate_pattern_match(seq_one, seq_two, int(mismatches))))


    # Approximate Pattern count
    seq_one = cfg.approximate_pattern_count_sample.pattern
    seq_two = cfg.approximate_pattern_count_sample.sequence
    mismatches = cfg.approximate_pattern_count_sample.mismatches
    print(approximate_pattern_count(seq_one, seq_two, int(mismatches)))


    for file in folder_txt_files(cfg.folders.approximate_pattern_count_sample):
        print(file)
        text = extract_seq_from_txt(file).split('\n')
        print(approximate_pattern_count(text[0], text[1], int(text[2])))


    print(immediate_neighbors("ATGC"))
    """
    # most frequent words with mismatches
    for file in folder_txt_files(cfg.folders.frequent_words_with_mismatches):
        print(file)
        text = extract_seq_from_txt(file).split('\n')
        # turn strings into integers and unpack values
        kmer_length , mismatches = [int(x) for x in text[1].split(" ")]
        frequent_words_with_mismatches(text[0], kmer_length, mismatches)
    """


if __name__ == "__main__":
    main()

    import doctest
    doctest.testmod()
