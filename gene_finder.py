# -*- coding: utf-8 -*-
"""
Created on Mon Sep  11 2017
Gene Finder Mini Project

@author: Shreya Rangarajan
"""

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
    """
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'G':
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
    """
    reverse_complement = ''
    for i in dna[::-1]:
        reverse_complement += get_complement(i)

    return reverse_complement

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
    """
    stop_codons = ['TAA', 'TAG', 'TGA']
    for i in range(0,len(dna)-1,3):
        seg = dna[i:i+3]
        if seg in stop_codons:
            dna = dna[0:i]
            return dna
    return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    new_ORFs_oneframe_list = []
    start_codon = 'ATG'
    i = 0
    while i < len(dna):
        seg = dna[i:i+3]
        if seg == start_codon:
            new_ORFs_oneframe_list.append(rest_of_ORF(dna[i:]))
            i += len(rest_of_ORF(dna[i:]))
        i += 3

    return new_ORFs_oneframe_list

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    new_ORFs_list_all = []
    for i in range(3):
        new_ORFs_list_all = new_ORFs_list_all + find_all_ORFs_oneframe(dna[i:])

    return new_ORFs_list_all

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse_complement_DNA = get_reverse_complement(dna)
    two_strands_list = find_all_ORFs(dna) + find_all_ORFs(reverse_complement_DNA)

    return two_strands_list

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest_ORF_string = max(find_all_ORFs_both_strands(dna))

    return longest_ORF_string


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest_ORF_overnumtrials = []
    for _ in range(num_trials):
        shuffle_DNA = shuffle_string(dna)
        longest_ORF_overnumtrials.append(len(longest_ORF(shuffle_DNA)))
    max_of_longest_orf = max(longest_ORF_overnumtrials)

    return max_of_longest_orf

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
    """
    amino_acids_list = ''
    index = 0

    if len(dna)%3 != 0:
        dna = dna[:-(len(dna)%3)]

    for index in range(0,len(dna),3):
        codon = dna[index:index+3]
        amino_acid_letter = aa_table[codon]
        amino_acids_list = amino_acids_list + amino_acid_letter
    return amino_acids_list

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna,1500)
    all_ORFs_bothstrands = []
    aa_sequence = []

    all_ORFs = find_all_ORFs_both_strands(dna)
    for orfs in all_ORFs:
        if len(orfs) > threshold:
            aa_sequence.append(coding_strand_to_AA(orfs))
    return aa_sequence

if __name__ == "__main__":
    import doctest
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
    #doctest.testmod()
    # doctest.run_docstring_examples(coding_strand_to_AA, globals())
