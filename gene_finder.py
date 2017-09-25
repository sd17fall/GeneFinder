# -*- coding: utf-8 -*-
"""
Finds amino acid sequences coded by DNA input

@author: Vivien Chen

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")
import doctest
from pickle import dump, load


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

        I added two more doctests because each complementary nucleotide is in
        its own if/else if branch. If one branch doesn't work and there is no
        doctest for that branch, the mistake will not be caught. In fact, I
        found out that I spelled nucleotide wrong in the G branch.

    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T':
        return 'A'

#doctest.run_docstring_examples(get_complement, globals(), verbose=True)


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        I did not add any more doctests because each string contains all the
        nucleotides, so two doctests with random strings is sufficient.

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reversed_dna = dna[::-1]
    comp_list = []

    for nuc in range(0, len(reversed_dna)):
        comp_nuc = get_complement(reversed_dna[nuc])
        comp_list.append(comp_nuc)

    comp_str = ''.join(comp_list)

    return comp_str

#doctest.run_docstring_examples(get_reverse_complement, globals(), verbose=True)


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        I added two more doctests: one in which the frame stop codon is TAA and
        another in which there is no frame stop codon. This way, all the
        branches are tested for error.

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGGAATAAGTA")
    'ATGGAA'
    >>> rest_of_ORF("ATGAGAGA")
    'ATGAGAGA'
    """
    stop_codons = ['TAA', 'TAG', 'TGA']

    for i in range(3, len(dna)-2, 3):
        if dna[i:i+3] in stop_codons:
            return dna[:i]

    return dna

#doctest.run_docstring_examples(rest_of_ORF, globals(), verbose=True)


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        I added another doctest that includes three nested ATGs in the default
        frame of the sequence in order to make sure that the function does not
        return nested ORFs in that frame.

    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe('ATGATGATGTAAC')
    ['ATGATGATG']
    """
    all_ORFs = []
    stopped = True
    stop_codons = ['TAA', 'TAG', 'TGA']

    for i in range(0, len(dna)-2, 3):
        if dna[i:i+3] == 'ATG' and stopped:
            all_ORFs.append(rest_of_ORF(dna[i:]))
            stopped = False
        elif dna[i:i+3] in stop_codons:
            stopped = True

    return all_ORFs

#doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=True)


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        I did not add any more doctests because this function depends on the
        previous functions, so assuming the doctests of those were fruitful,
        one doctest is sufficient for this function. The doctest includes one
        ORF in each frame, so the function works if the doctest yields true.

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    frame1 = find_all_ORFs_oneframe(dna)
    frame2 = find_all_ORFs_oneframe(dna[1:])
    frame3 = find_all_ORFs_oneframe(dna[2:])

    return frame1 + frame2 + frame3

#doctest.run_docstring_examples(find_all_ORFs, globals(), verbose=True)


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        I did not add any more doctests because this function depends on the
        previous functions, so assuming the doctests of those were fruitful,
        one doctest is sufficient for this function. The doctest includes one
        ORF in each strand, so the function works if the doctest yields true.

    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    strand1 = find_all_ORFs(dna)
    strand2 = find_all_ORFs(get_reverse_complement(dna))

    return strand1 + strand2

#doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose=True)


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_ORFs = find_all_ORFs_both_strands(dna)
    humongest = 0

    for i in range(0, len(all_ORFs)-1):
        if len(all_ORFs[i]) < len(all_ORFs[i+1]):
            humongest = all_ORFs[i+1]

    return humongest

#doctest.run_docstring_examples(longest_ORF, globals(), verbose=True)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    dna_shuffles = []

    for i in range(0, num_trials):
        dna_shuffles.append(shuffle_string(dna))

    for i in range(0, len(dna_shuffles)-1):
        if longest_ORF(dna_shuffles[i]) < longest_ORF(dna_shuffles[i+1]):
            humongest = longest_ORF(dna_shuffles[i+1])

    return len(humongest)

#print longest_ORF_noncoding('ATGCGAATGTAGCATCAAA', 100)


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
    amino_acids = []

    for i in range(0, len(dna)-2, 3):
        amino_acid = aa_table[dna[i:i+3]]
        amino_acids.append(amino_acid)

    amino_acids_str = ''.join(amino_acids)

    return amino_acids_str

#doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    all_ORFs = find_all_ORFs_both_strands(dna)
    amino_acids = []

    for i in range(0, len(all_ORFs)):
        if len(all_ORFs[i]) > threshold:
            amino_acid = coding_strand_to_AA(all_ORFs[i])
            amino_acids.append(amino_acid)

    return amino_acids


if __name__ == "__main__":
    #import doctest
    #doctest.testmod()

    genes = gene_finder(dna)

    #write to genes.txt file
    file_ = open('genes.txt', 'wb')
    dump(genes, file_)

    #read from genes.txt file
    file_ = open('genes.txt', 'rb+')
    print(load(file_))
