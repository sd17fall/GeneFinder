# -*- coding: utf-8 -*-
"""
First project for Olin Software Design Fall 2017

@author: Emma Westerhoff

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
    nucleotide_inputs = ['A', 'T', 'C', 'G']
    nucleotide_complements = ['T', 'A', 'G', 'C']
    i = 0
    complement = 'x' #lets the user know the complement was incorrectly computed
    while i < len(nucleotide_inputs):
        if nucleotide_inputs[i] == nucleotide:
            complement = nucleotide_complements[i]
        i += 1
    return complement

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
    reverse = ''
    i = 0
    length = len(dna)

    while i < length:
        letter = dna[length - 1 -i] #moves backwards along the string
        pair = get_complement(letter) #finds complement
        reverse = reverse + pair
        i += 1
    return reverse

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
    >>> rest_of_ORF("ATTTCGGGT")
    'ATTTCGGGT'
    """
    stop_codons = ['TAG', 'TGA', 'TAA']
    codons = []
    n = 3

    for i in range(0, len(dna), n):
        codons.append(dna[i:i+n])

    for c in range(0, len(codons)):
        for s in range(0, len(stop_codons)):
            if codons[c] == stop_codons[s]:
                codons = codons[:c]
                return_string = ''.join(codons)
                return return_string

    return_string = ''.join(codons)
    return return_string

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
    codons = []
    n = 3
    ORFS = []
    c = 0

    for i in range(0, len(dna), n):
        codons.append(dna[i:i+n])

    while c in range(0, len(codons)):
        if codons[c] == start_codon:
                dna_sequence = rest_of_ORF(''.join(codons[c:]))
                ORFS.append(dna_sequence)
                c += len(dna_sequence) #skips over the rest of the sequence
        c += 1 #if I'm missing a permutation, this might be a problem.

    return ORFS


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
    return_list = []

    for i in range (0,3):
        #cases are coming through that occur in the same frame
        orfs = find_all_ORFs_oneframe(dna[i:])
        for o in orfs:
            result = ''.join(o)
            if result != '': #if there are no permutations in a run through
                return_list.append(result)
    return return_list

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    return_list = []

    return_list.append(''.join(find_all_ORFs(dna)))
    return_list.append(''.join(find_all_ORFs(get_reverse_complement(dna))))
    return return_list


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orfs = find_all_ORFs_both_strands(dna)

    longest_size=len(max(orfs,key=len))

    for o in orfs:
        if(longest_size==len(o)):
            longest = o

    return longest

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    lengths = []
    for rand in range(0, num_trials):
        new_sequence = shuffle_string(dna)
        leng = len(longest_ORF(new_sequence))
        lengths.append(leng)

    maximum = max(lengths)
    return maximum


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
    n = 3
    codons = []
    acids = ''

    for i in range(0, len(dna), n):
        codons.append(dna[i:i+n])

    if len(codons[-1]) < 3:
        codons.pop(-1)

    for c in codons:
        amino = aa_table[c]
        acids += ''.join(amino)

    return acids


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 15)
    #change this to 1500 or so later
    dna_orfs = find_all_ORFs_both_strands(dna)
    amino_sequences = []
    longs = []

    for snip in dna_orfs:
        if len(snip) > threshold:
            longs.append(snip)

    for string in dna_orfs:
        amino_sequences.append(coding_strand_to_AA(''.join(string)))

    #print(amino_sequences)
    #print(len(amino_sequences))
    return amino_sequences

if __name__ == "__main__":
    import doctest
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
    dna = load_seq("./data/X73525.fa")
    gene_finder(dna)
