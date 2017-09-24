# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

dna = load_seq("./data/X73525.fa")


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
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'C':
        return 'G'
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
    complement_list= []                          #creating empty complement list
    for letter in dna:                           #creating for loop to check every letter in the string
        complement = get_complement(letter)      #getting complements of string
        complement_list.append(complement)       #complement list filld with complements(adding complement list to empty list)
        # reverse the complement list
    complement_list.reverse()                    #reversing the complement list
    return ''.join(complement_list)                     #finsl step is to make the complement list into one string


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
    stop_codon_list = ['TAG', 'TAA', 'TGA']
    orf = ''
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if codon in stop_codon_list:
            return orf
        else:
            orf = orf + codon
    return dna[:]


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
    start_codon = 'ATG'
    orf_list = []
    i=0
    while i < len(dna):
        codon = dna[i:i+3]
        if codon == start_codon:
            orf = rest_of_ORF(dna[i:])
            orf_list.append(orf)
            i = i + len(orf)
        else:
            i = i + 3
    return orf_list


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

    a = find_all_ORFs_oneframe(dna)
    b = find_all_ORFs_oneframe(dna[1:])
    c = find_all_ORFs_oneframe(dna[2:])
    return a + b + c


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    d = find_all_ORFs(dna)
    e = find_all_ORFs(get_reverse_complement(dna))
    return d + e


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    max_ORF = ''
    max_ORF = max(find_all_ORFs_both_strands(dna), key=len)
    return max_ORF


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    shuffeled = ''
    longest_shuffeled_ORF = ''
    shuffeled_list = []
    i= 0
    while i < num_trials:
        shuffeled = shuffle_string(dna)
        longest_shuffeled_ORF = longest_ORF(shuffeled)
        shuffeled_list.append(longest_shuffeled_ORF)
        i = i + 1
    return len(max(shuffeled_list, key=len))


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
    amino_acid_list = []
    i = 0
    while i < len(dna):
        codon = dna[i:i+3]
        if len(dna[i:]) < 3:
            return ''.join(amino_acid_list)
        amino_acid = aa_table[codon]
        amino_acid_list.append(amino_acid)
        i = i + 3
    return ''.join(amino_acid_list)



def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        Next, find all open reading frames on both strands, and then return a list
        containing the amino acid sequence encoded by any open reading frames that
        are longer than the threshold computed above using longest_ORF_noncoding.


    """
    amino_acids_longer_than_threshold_string = ''
    amino_acids_longer_than_threshold = []
    threshold = longest_ORF_noncoding(dna, 1500)
    all_ORFs_both_stands = find_all_ORFs_both_strands(dna)
    for ORF in all_ORFs_both_stands:
        if(len(ORF) > threshold):
            amino_acids_longer_than_threshold_string = coding_strand_to_AA(ORF)
            amino_acids_longer_than_threshold.append(amino_acids_longer_than_threshold_string)
    return  amino_acids_longer_than_threshold

if __name__ == "__main__":
    import doctest

    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
    #doctest.testmod(verbose = True)
    doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
    #longest_ORF_noncoding(dna,500)
