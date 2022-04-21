"""
Jim DelProposto
Bioinformatics Specialization
Week 1
Finding Hidden Messages in DNA (Bioinformatics I)
"""


import urllib.request
import os



def extract_seq_from_txt(file_path: os.path) -> str:

    """
    extracts a string sequence from a text file
    param: file_path: path where the file is located
    return: sequence: sequence extracted from the text file

    >>> len(extract_seq_from_txt(os.path.normpath(".\data\Vibrio_cholerae.txt")))
    1108250
    >>> extract_seq_from_txt(os.path.normpath(".\data\bad_loc.txt"))
    -1
    """


    if os.path.exists(file_path):
        text_file = open(file=file_path, mode="r", encoding='utf-8')
        sequence = text_file.read()
    else:
        sequence = -1

    return sequence

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
            text_file = open(file_path, "w")
            text_file.write(sequence)
            text_file.close
        except:
            print("unable to open URL")
            return -1

    return os.path.normpath(file_path)


def main():
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
if __name__ == "__main__":
    main()
    import doctest
    doctest.testmod



