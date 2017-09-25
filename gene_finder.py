# -*- coding: utf-8 -*-
"""
ENGR2510 - Software Design - Fall 2017
This program analyzes a DNA sequence and outputs snippets of DNA that are likely
to be protein-coding genes.
@author: Victor Bianchi
"""
#FINALLLLLLLLLLLLLLL

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    # TODO: implement this
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    else:
        return 'C'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("ATCG")
    'CGAT'
    """
    # TODO: implement this
    reversed_dna = dna[::-1]
    result = ' '
    for letter in reversed_dna:
        result = result + get_complement(letter)
        return result

    def divide_to_codons(dna):
        """Takes a DNA sequence and outputs a list of string triplets(codons) that makes up the sequence
           Last element might be incomplete codon with less then three letters
        >>> divide_to_codons("ATGTGAA")
        ['ATG', 'TGA', 'A']
        >>> divide_to_codons("ATGTGA")
        ['ATG', 'TGA']
        >>> divide_to_codons("ATGTGAAA")
        ['ATG', 'TGA', 'AA']
        """
        index = 0
        result = []
        while index < len(dna):
            result.append(dna[index:index+3])
            index = index + 3
        return result

    def rest_of_ORF(dna):
        """ Takes a DNA sequence that is assumed to begin with a start
            codon and returns the sequence up to but not including the
            first in frame stop codon.  If there is no in frame stop codon,
            returns the whole string.
            dna: a DNA sequence
            returns: the open reading frame represented as a string
        >>> rest_of_ORF("ATGTGAA")
        'ATG'
        >>> rest_of_ORF("ATGAGATAGG")
        'ATGAGA'
        >>> rest_of_ORF("ATG")
        'ATG'
        >>> rest_of_ORF("AT")
        'AT'
        >>> rest_of_ORF("ATGASDASDWASDWADASDSAD")
        'ATGASDASDWASDWADASDSAD'
        >>> rest_of_ORF("ATGTGTTAAATGAAAAAATAGAA")
        'ATGTGT'
        """
        stop_codons = ['TAG', 'TAA, TGA']
        #list of codons from which the dna is composed of
        codons = divide_to_codons(dna)
        result = ""
        index = 0
        while index + 1 < len(codons):
            #If next codons isn't a stop codon, add it to string and iterate
            if codons[index + 1] not in stop_codons:
                result = result + codons[index]
                index = index + 1
            else:
                #Add codon before stop codon
                result = result + codons[index]
                return result
        return dna

    def find_all_ORFs_oneframe(dna):
        """
        Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
        >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
        ['ATGCATGAATGTAGA', 'ATGTGCCC']
        >>> find_all_ORFs_oneframe("ATGTGAA")
        ['ATG']
        >>> find_all_ORFs_oneframe('ASDASDAWSDSD')
        []
        >>> find_all_ORFs_oneframe('TATATGCATGAATGTAGATAGATGTGCTAAATAATAATGTTTTAAATT')
        ['ATGCATGAATGTAGA', 'ATGTGC', 'ATGTTT']
        """
        index = 0
        orf_list = []
        while index < len(dna):
            if dna[index:index+3] == 'ATG':
                #appended ORF
                orf = rest_of_ORF(dna[index:])
                orf_list.append(orf)
                index = index + len(orf)
            else:
                index = index + 3
        return orf_list

    def find_all_ORFs(dna):
        """
        Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
        This unit testing would be enough because there isn't any special exceptions that needs to be tested. Also, this case tests this function's
        ability to grab orf from three different possible reading frames.
        >>> find_all_ORFs("ATGCATGAATGTAG")
        ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
        """
        #orf list from all frames
        orf_list = []
        #zero offset frame
        orf_list = orf_list + find_all_ORFs_oneframe(dna)
        #first offset frame
        orf_list = orf_list + find_all_ORFs_oneframe(dna[1:])
        #second offset frame
        orf_list = orf_list + find_all_ORFs_oneframe(dna[2:])
        return orf_list

    def find_all_ORFs_both_strands(dna):
        """
        Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
        >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
        ['ATGCGAATG', 'ATGCTACATTCGCAT']
        """
        reverse = get_reverse_complement(dna)
        #finds orfs in both direction
        orf_list = find_all_ORFs(dna) + find_all_ORFs(reverse_complement)
        return orf_list

    def longest_ORF(dna):
        """
        Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        >>> longest_ORF("ATGCGAATGTAGCATCAAA")
        'ATGCTACATTCGCAT'
        """
        longest_length = 0
        orfs = find_all_ORFs_both_strands(dna)
        for orf in orfs:
            if len(orf) > longest_length:
                longest_orf = orf
                longest_length = len(orf)
            return longest_orf

    def longest_ORF_noncoding(dna, num_trials):
        """
        Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        """
        x = 0
        longest = 0
        while x < num_trials:
            shuffled_dna = shuffle_string(dna)
            longest_orf_length = len(longest_ORF(shuffled_dna))
            if longest_orf_length > longest:
                longest = longest_orf_length
            x = x + 1
        return longest

        def coding_strand_to_AA(dna):
            """ Computes the Protein encoded by a sequence of DNA.  This function
                does not check for start and stop codons (it assumes that the input
                DNA sequence represents an protein coding region).
                dna: a DNA sequence represented as a string
                returns: a string containing the sequence of amino acids encoded by the
                         the input DNA fragment
                >>> coding_strand_to_AA("ATGCGA")
                'MR'
                >>> coding_strand_to_AA("ATGCCCGCTTT")
                'MPA'
                >>> coding_strand_to_AA("TTTATCATGTTAGTTA")
                'FIMLV'
            """
            codons = divide_to_codons(dna)
            amino_acid = ''
            for codon in codons:
                if len(codon) == 3:
                    amino_acid = amino_acid + aa_table[codon]
            return amino_acid

        def gene_finder(dna):
            """ Returns the amino acid sequences that are likely coded by the specified dna
                dna: a DNA sequence
                returns: a list of all amino acid sequences coded by the sequence dna.
            """
            threshold = longest_ORF_noncoding(dna, 1500)
            all_orfs = find_all_ORFs_both_strands(dna)
            amnio_acids = []
            for orf in all_orfs:
                if len(orf) > threshold :
                    amino_acids.append(coding_strang_to_AA(orf))
            return amino_acids

    if __name__ == "__main__":
        import doctest
        doctest.testmod(verbose = True)
        doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose = True)
        dna_seq = load_seq('data/X73525.fa')
        print (gene_finder(dna_seq))
