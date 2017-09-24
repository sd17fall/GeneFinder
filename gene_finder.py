# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Harris Davidson

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
    if nucleotide == 'A': return 'T'
    elif nucleotide == 'T': return 'A'
    elif nucleotide == 'C': return 'G'
    elif nucleotide == 'G': return 'C'


# print(get_complement('A'))
# print(get_complement('T'))
# print(get_complement('C'))
#print(get_complement('G'))

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
    compliment = ""
    for f in dna:
        compliment = compliment + get_complement(f)
    return compliment[::-1]

#print(get_reverse_complement("ATGCCCGCTTT"))

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
    i = 0
    while i < len(dna):
        codon = dna[i:i+3]
        if (codon == 'TAG') or (codon == 'TAA') or (codon == 'TGA'):
            return dna[0:i]
        i = i+3
    return dna

#print(rest_of_ORF("ATGAGATAGG"))
    # pos = [] #indeices of stop codons
    # b = 0
    # while b>-1:
    #     b = dna.find("TAG",b+1)
    #     pos = pos + [b]
    # b=0
    # while b>-1:
    #     b = dna.find("TAA",b+1)
    #     pos = pos + [b]
    # b=0
    # while b>-1:
    #     b = dna.find("TGA",b+1)
    #     pos = pos + [b]
    #
    # stops = []
    #
    # for x in pos:
    #     if 0==x%3:
    #         stops.append(x)
    # stops.sort
    # if len(stops) == 0:
    #     return dna
    # else:
    #     return dna[0:stops[-1]]

#print(rest_of_ORF('ATGAGATAGG'))

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
    all_ORFs = []
    #print(all_ORFs)

    #for i in range(0,len(dna),3):
    i = 0
    while i < len(dna)-3:
        if dna[i:i+3] == 'ATG':
            #print(rest_of_ORF(dna[i:]))
            all_ORFs.append(rest_of_ORF(dna[i:]))
            #print(all_ORFs[-1])
            #print(len(all_ORFs[-1]))
            i = i + len(all_ORFs[-1])
            #print(i)
        else:
            i = i + 3
    return all_ORFs
#print(find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC"))

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
    return find_all_ORFs_oneframe(dna)+find_all_ORFs_oneframe(dna[1:])+find_all_ORFs_oneframe(dna[2:])
#print(find_all_ORFs('ATGCATGAATGTAG'))

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
#print(find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA"))

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_ORFs = find_all_ORFs_both_strands(dna)
    length = []
    for f in all_ORFs:
        length.append(len(f))
    longest = length.index(max(length))
    return all_ORFs[longest]
    #print(all_ORFs[longest])


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    max_length = 0
    length = []
    for i in range(1,num_trials):
        sdna = shuffle_string(dna)
        if(max_length < len(longest_ORF(sdna))):
            max_length = len(longest_ORF(sdna))
    return max_length

#print(longest_ORF_noncoding('ATGCGAATGTAGCATCAAA',3))

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
    acid_sequence = ''
    for i in range(0,len(dna)-2,3):
        #print(type(i))
        #codon = dna[i:i+3]

        #print(codon)
        acid_sequence =  acid_sequence + aa_table[dna[i:i+3]]
    #print(acid_sequence)
    return acid_sequence

#print(coding_strand_to_AA("ATGCCCGCTTT"))


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    minlength = longest_ORF_noncoding(dna,1500)
    #minlength = 6000
    print(minlength)
    all_ORFs = find_all_ORFs_both_strands(dna)
    #print(all_ORFs)
    all_coding_ORFs=[]
    for i in range(0,len(all_ORFs)-1):
        if len(all_ORFs[i]) > minlength:
            #print(all_ORFs[i])
            all_coding_ORFs.append(all_ORFs[i])
        #else:
            #print('too short')
    #print(all_coding_ORFs)
    amino_sequences=[]
    for i in range(0,len(all_coding_ORFs)):
        #print(coding_strand_to_AA(all_coding_ORFs[i]))
        #print(all_coding_ORFs[i])
        amino_sequences.append(coding_strand_to_AA(all_coding_ORFs[i]))

    return amino_sequences

#print(gene_finder("ATGCCCGCTTT"))

from load import load_seq
dna = load_seq("./data/X73525.fa")

print(gene_finder(dna))
#gene_finder(dna)
#
# if __name__ == "__main__":
#     import doctest
#     doctest.testmod(verbose=True)
