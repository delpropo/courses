"""
Jim DelProposto
Bioinformatics Specialization
Week 1
Finding Hidden Messages in DNA (Bioinformatics I)
"""



import urllib.request
import os
from collections import OrderedDict

import numpy as np

COMPLEMENT_DNA_BASES = {"A": "T", "T": "A", "G": "C", "C": "G"}
DNA_BASES = ("A", "T", "G", "C")

def extract_seq_from_txt(file_path: os.path) -> str:
    """
    Extract string from file.  This may include multiple strings.
    removes any new line from end of string.
    extracts a string sequence from a text file
    param: file_path: path where the file is located
    return: sequence: string extracted from the text file

    >>> len(extract_seq_from_txt(os.path.normpath(".\data\Vibrio_cholerae.txt")))
    1108250
    >>> extract_seq_from_txt(os.path.normpath(".\data\bad_loc.txt"))
    -1
    """


    if os.path.exists(file_path):
        with open(file=file_path, mode="r", encoding='utf-8') as text_file:
            return text_file.read().rstrip('\n')
    else:
        return -1



def pattern_count(sequence: str, pattern: str) -> int:
    """
    search a sequence for a pattern

    return the number of times that pattern appears
    param: sequence:  string with the sequence being searched
    param: pattern:  string being searched for in the pattern
    return: count of the number of patterns found

    >>> pattern_count("GCGCG", "GCG")
    2
    """

    count = 0
    pattern_len = len(pattern)
    for i in range(len(sequence)):
        if pattern == sequence[i:(i + pattern_len)]:
            count += 1
    return count

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

def pattern_locations(sequence: str, pattern: str) -> list:
    """
    searches a sequences for the starting bp locations
    location start with bp as 0, not 1
    param: sequence:  string with the sequence being searched
    param: pattern:  the sequences being searched
    return: a list of the bp locations of the start of the sequence
    >>> " ".join([str(j) for j in pattern_locations("GATATATGCATATACTT", "ATAT")])
    '1 3 9'
    """
    locations = []
    pattern_len = len(pattern)
    for i in range(len(sequence)):
        if pattern == sequence[i:(i + pattern_len)]:
            locations.append(i)
    return locations

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

def freq_set_search(sequence: str, seq_set: set) -> dict:
    """
    search for sequences in a set and count the times they occur

    param: sequence:  string with the sequence being searched
    param: pattern:  string being searched for in the pattern
    return: count of the number of patterns found


    """
    freq_dict = {}
    for i in seq_set:
        if i in freq_dict:
            raise ValueError(f"{i} is contained in freq dict and has not been counted")
        freq_dict[i] = pattern_count(sequence, i)

    return freq_dict

def max_freq_words(sequence: str, size: int, all_values = False ) -> OrderedDict:
    """
    search for sequences of specific sizes and count the times they occur
    Created for second code challenge
    param: sequence:  string with the sequence being searched
    param: pattern:  string being searched for in the pattern
    param:  all_values can be set to True is an ordered dict with all values is required
    return: the patterns and number of times they occur

    >>> " ".join(max_freq_words("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4).keys())
    'GCAT CATG'
    >>> " ".join(max_freq_words("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, all_values = True).keys())
    'ACGT CGTT GTTG TTGC TGCA GCAT CATG ATGT TGTC GTCG TCGC CGCA ATGA TGAT GATG ATGC TGAG GAGA AGAG GAGC AGCT'
    """
    seq_dict = OrderedDict()
    max_words = OrderedDict()
    for i in seq_size_generator(sequence, size):
        assert len(i) == size
        if i in seq_dict.keys():
            seq_dict[i] += 1
        elif i is not None:
            seq_dict[i] = 1
    max_freq = max(seq_dict.values())

    # extract only values which match the highest frequency
    if all_values:
        return seq_dict

    for (key, value) in seq_dict.items():
        if value == max_freq:
            max_words[key] = value

    return max_words



def dl_save_txt(url_str: str, folder=r".\data") -> os.path:
    """
    download a text file
    save it to the data folder if the file does not exist
    must be run from current working directory

    VIBRIO_CHOLERAE_LENGTH = 1108250

    vb_genome_url = "http://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt"
    file_path = ".\data\Vibrio_cholerae.txt"
    """

    txt_name = url_str.split("/")[-1]
    file_path = os.path.join(os.path.normpath(folder), txt_name)

    if not os.path.exists(file_path):
        try:
            data = urllib.request.urlopen(url_str)
            sequence = str(data.read(), 'utf-8')
             # save downloaded file
            with open(file_path, "w") as text_file:
                text_file.write(sequence)
        except ValueError:
            print("unable to open URL")
            return -1

    return os.path.normpath(file_path)

def reverse_complement(sequence: str) -> str:
    """
    returns the reverse complement of a sequence
    example of  ATGC -> GCAT

    param: sequence: The DNA string being examined
    return: the reverse complement of the sequence
    >>> reverse_complement("ATGC")
    'GCAT'
    >>> reverse_complement("AAAAAAAAAA")
    'TTTTTTTTTT'
    >>> reverse_complement("ATGTAATAG")
    'CTATTACAT'
    """
    return "".join([COMPLEMENT_DNA_BASES[i] for i in sequence[::-1] ])

def search_clumps(sequence: str, size: int, range_clumps: int, min_clumps: int ) -> set():
    """
    Search for groupings of sequences that fall within a specified range

    param: sequence: The DNA string being examined
    param: size: The size of the k-mer being looked for
    param: range_clumps: The sequence range being looked at in bp
    param: min_clumps: The minimum number of sequence clumps in the specified range
    return: the reverse complement of the sequence
    >>> sorted(search_clumps("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4))
    ['CGACA', 'GAAGA']
    """
    # get start values for everything
    # if the start values in a specific range are > a certain value then success

    # creates an ordered dict with all sequences and the number of times they occur
    frequency_dict = max_freq_words(sequence, size, all_values=True)

    # First round of removal.  remove any sequnces that show up less than min_occurances
    #search through all max_freq_words dict if value is greater or equal to min_occurances



    pattern_list = set()
    for (key, value) in frequency_dict.items():
        if value >= min_clumps:
            #start_list = pattern_locations(sequence, key)
            if max_values_in_range(pattern_locations(sequence, key), range_clumps) >= min_clumps:
                pattern_list.add(key)

    return pattern_list




def max_values_in_range(value_list: list, bp_range: int) -> int:
    """
    Determines the maximum number of values in a specified range of a list
    [0, 1, 2, 5, 10] with a range of 2 would have 2 values

    param: value_list: list of values which should be in numerical order
    param: bp_range the largest difference between two values in the list allowed
    return:  max_value:  the largest number of values in a slice where the starting and final value are less than or equal to range_num
    >>> values_in_range([0, 10, 50, 85, 86, 91], 5)
    2
    >>> max_values_in_range([0, 10, 50, 84, 85, 86, 91], 6)
    3
    >>> max_values_in_range([0, 1, 2, 3, 4, 5, 6], 6)
    6
    >>> max_values_in_range([1, 2, 3, 4, 5, 6], 6)
    6
    >>> max_values_in_range([0, 499, 500], 6)
    2
    """
    #sort(values_in_range)
    num_values = len(value_list)
    # create a default of 0 for the maximum number of starts is a given bp_range
    max_value = 0
    for i in range(num_values):
        for j in range(i+1, num_values):
            assert j > i and value_list[j] > value_list[i]
            if bp_range > (value_list[j] - value_list[i]) and max_value < j-i + 1:
                max_value = j-i + 1

    return max_value



def main() -> None:
    # vibrio cholerae genome is 1108250 nucleotides
    VIBRIO_CHOLERAE_LENGTH = 1108250
    vb_genome_url = "http://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt"
    file_path = ".\data\Vibrio_cholerae.txt"
    sequence = extract_seq_from_txt(dl_save_txt(vb_genome_url))
    assert (len(sequence)) == VIBRIO_CHOLERAE_LENGTH

    test_file_dir = os.path.normpath(".\data\PatternCount\inputs")
    file_list = [os.path.join(test_file_dir, file) for file in os.listdir(test_file_dir) if file.endswith(".txt")]
    file_list.sort()
    for i in file_list:
        print(i)
        sequences = extract_seq_from_txt(i).split("\n")
        print(pattern_count(sequences[0], sequences[1]))

    # test case
    values_freq = max_freq_words("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)
    print(" ".join(values_freq.keys()))

    # if sorting is required.  Manually pulled information and pasted it due to no example txt file starting out
    values_freq = max_freq_words("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)
    print(" ".join(sorted(values_freq.keys())))


    # reverse complement
    reverse_complement("ATCGCAGAC")


    test_file_rc_dir = os.path.normpath(".\data\ReverseComplement\inputs")
    file_list = [os.path.join(test_file_rc_dir, file) for file in os.listdir(test_file_rc_dir) if file.endswith(".txt")]
    file_list.sort()
    for i in file_list:
        sequences = extract_seq_from_txt(i).split("\n")
        print("reverse complement sequence")
        print(reverse_complement(sequences[0]))

    # locate the starting position of a pattern in a sequence.  Space separated

    test_file_start_dir = os.path.normpath(".\data\PatternMatching_locations\inputs")
    file_list = [os.path.join(test_file_start_dir, file) for file in os.listdir(test_file_start_dir) if file.endswith(".txt")]
    file_list.sort()
    for i in file_list:
        sequences = extract_seq_from_txt(i).split("\n")
        print(i)
        print(" ".join([str(j) for j in pattern_locations(sequences[1], sequences[0])]))

    # process the vibrio cholorea genome through pattern_locations
    cholera_file = os.path.normpath(".\data\PatternMatching_locations\inputs\Vibrio_cholerae.txt")
    print(cholera_file)
    sequences = extract_seq_from_txt(cholera_file).split("\n")
    locations = " ".join([str(j) for j in pattern_locations(sequences[0], "CTTGATCAT")])
    print(locations)

    # clump search
    test_file_start_dir = os.path.normpath(".\data\clump_search")
    file_list = [os.path.join(test_file_start_dir, file) for file in os.listdir(test_file_start_dir) if file.endswith(".txt")]
    file_list.sort()
    for i in file_list:
        sequences = extract_seq_from_txt(i).split("\n")
        #print(search_clumps(sequences, *[int(i) for i in sequences[1].split(" ")]))
        # value ended up being an empty set



if __name__ == "__main__":
    main()
    import doctest
    doctest.testmod
