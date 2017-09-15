# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: JOHN WEN

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
    """ Checks for each letter and returns the complement
        (might not be the most effiicent way) and returns
        'This is not a nucleotide' if any other letters are implented
        Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'This is not a nucleotide'


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
    reverse_dna = ''
     # base string to add complements
    index = len(dna) - 1
     # counter for while loop
    while index > 0:
        reverse_dna = reverse_dna + get_complement(dna[index])
        index = index - 1
    #  while loop that processes each letter of the string and applies the
    #  previous function get_complement() to add on to the newly made string
    #  reverse_dna
    reverse_dna = reverse_dna + get_complement(dna[index])
    #adds the complement of the last element of the string
    return reverse_dna


#This one code took me an hour to write. jesus.
def rest_of_ORF(dna):
    """The rest of the function checks every three codons and tries to match them
    to the stop codon. If it isn't a stop codon then it will add the letter to the string

    Takes a DNA sequence that is assumed to begin with a start
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
    index = 0
    newstr = ''
    stop_codon = ['TAG','TAA','TGA']
    if ('TAG' not in dna) and ('TAA' not in dna) and ('TGA' not in dna):
        return dna
    else:
        while index < len(dna):
            if (dna[index:index+3] not in stop_codon):
                newstr = newstr + dna[index:index+3]
                index = index + 3
            else:
                return newstr
        return newstr

# V1
# def find_all_ORFs_oneframe(dna):
#   newstr1 = []
#   newstr1.append(rest_of_ORF(dna))
#   stop_codon = ['TAG','TAA','TGA']
#   index = 0
#   #newstr1 adds the first chunk to the list
#   while dna[index:index+3] not in stop_codon:
#     index = index + 3
#   newstr1.append(dna[index+3:])
#   #checks to find
  # return newstr1


#V2 all_ORFS_oneframe
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
    >>> find_all_ORFs_oneframe("TTTATGCCCTAGATAATGTTTTAGATGCCCTAG")
    ['ATGCCC', 'ATGTTT', 'ATGCCC']
    >>> find_all_ORFs_oneframe("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATGTAG']
    """
    #print(find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC"))
    #print(find_all_ORFs_oneframe('AAAAAAAATAGAAA'))
    #print(find_all_ORFs_oneframe('AAAATTTTTTAATTTTTTATTATA'))
    finallistofORFS = []
    stop_codon = ['TAG','TAA','TGA']
    index = 0
    #newstr1 adds the first chunk to the list
    while index < len(dna):
        if dna[index:index+3] == 'ATG':
            index2 = index
            while index2 < len(dna):
                if dna[index2:index2+3] in stop_codon or index2 > len(dna) - 1:
                    finallistofORFS.append(dna[index:index2])
                index2 = index2 + 3
        index = index + 3
    return finallistofORFS





def find_all_ORFs(dna):
    """Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG']
    """
    finallistofORFs = []
    index = 0
    while index < 3:
        finallistofORFs.extend(find_all_ORFs_oneframe(dna[index:]))
        index = index + 1
    return finallistofORFs


#
def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    'TTTGATGCTACATTCGCAT'

    """
    finallistofORFs = []
    finallistofORFs.extend(find_all_ORFs(dna))
    finallistofORFs.extend(find_all_ORFs(get_reverse_complement(dna)))
    return finallistofORFs
#
#
# def longest_ORF(dna):
#     """ Finds the longest ORF on both strands of the specified DNA and returns it
#         as a string
#     >>> longest_ORF("ATGCGAATGTAGCATCAAA")
#     'ATGCTACATTCGCAT'
#     """
#     # TODO: implement this
#     pass
#
#
# def longest_ORF_noncoding(dna, num_trials):
#     """ Computes the maximum length of the longest ORF over num_trials shuffles
#         of the specfied DNA sequence
#
#         dna: a DNA sequence
#         num_trials: the number of random shuffles
#         returns: the maximum length longest ORF """
#     # TODO: implement this
#     pass
#
#
# def coding_strand_to_AA(dna):
#     """ Computes the Protein encoded by a sequence of DNA.  This function
#         does not check for start and stop codons (it assumes that the input
#         DNA sequence represents an protein coding region).
#
#         dna: a DNA sequence represented as a string
#         returns: a string containing the sequence of amino acids encoded by the
#                  the input DNA fragment
#
#         >>> coding_strand_to_AA("ATGCGA")
#         'MR'
#         >>> coding_strand_to_AA("ATGCCCGCTTT")
#         'MPA'
#     """
#     # TODO: implement this
#     pass
#
#
# def gene_finder(dna):
#     """ Returns the amino acid sequences that are likely coded by the specified dna
#
#         dna: a DNA sequence
#         returns: a list of all amino acid sequences coded by the sequence dna.
#     """
#     # TODO: implement this
#     pass

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
    # doctest.run_docstring_examples(find_all_ORFs_oneframe, globals())
    doctest.testmod()
