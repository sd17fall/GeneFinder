# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Felix Eberhardt

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


def get_complement(nucleotide): #week1
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
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'


def get_reverse_complement(dna): #week1
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # Step 1 get complement
    counter_1=0
    complement_dna = ''
    while counter_1 < len(dna):
        complement_dna = complement_dna + get_complement(dna[counter_1])
        counter_1 = counter_1 + 1
    # Step 2 reverse it
    reverse_complement = ''
    counter_2=len(dna)-1
    while counter_2 >= 0:
        reverse_complement = reverse_complement + complement_dna[counter_2]
        counter_2 = counter_2 - 1
    return reverse_complement

# Define Stop Codons
stop_codons = ['TAA', 'TAG', 'TGA']
start_codon = 'ATG'

def rest_of_ORF(dna): #week1
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
    i=0
    for i in range(0,len(dna), 3):
        if dna[i:i+3] in stop_codons:
            return dna[:i]
    return dna

def find_all_ORFs_oneframe(dna): #week1
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

    count=0
    all_ORFs_oneframe = [] #empty list
    while count < len(dna):
        if dna[count:count+3] == start_codon:
            orf = rest_of_ORF(dna[count:])
            all_ORFs_oneframe.append(orf)
            count = count + len(orf)
        else :
            count = count + 3
    return all_ORFs_oneframe

def find_all_ORFs(dna): #week1
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
    all_ORFs = find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna[1:]) + find_all_ORFs_oneframe(dna[2:])
    return all_ORFs

def find_all_ORFs_both_strands(dna): #week1
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
#   create empty lists
    all_ORFs_both_strands = []
#   search in them and add lists together
    all_ORFs_both_strands = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return all_ORFs_both_strands

    ##### week 2 ######

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest_ORF = max(find_all_ORFs_both_strands(dna), key=len)
    return longest_ORF

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    t=0
    ORF_noncoding = []
    while t < num_trials:        # loop it for num_trials times
        a = shuffle_string(dna) # use shuffling function to create new dna
        b = longest_ORF(a)      # put it in longest_ORF
        ORF_noncoding.append(b) # return len as integer to list
        t = t + 1
    longest_ORF_noncoding = max(ORF_noncoding, key=len)# Look for and return the longest ORF (same function as used above)
    max_length = len(longest_ORF_noncoding)
    return max_length



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

    strand_to_AA = '' # empty string to append
    for n in range(0, len(dna), 3): #go through dna string
        for i in range(0, len(codons)): #go through inner lists
            if dna[n:n+3] in codons[i]: #search for dna in each inn
                amino =  aa[i]   # convert each triplet into the linked letter
                strand_to_AA += amino # append it to the amino string
    return strand_to_AA #return the string of aminos


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500) # sets the treshold
    open_reading_frames = find_all_ORFs_both_strands(dna) # find all open reading frames on both strands
    i = 0
    aa_sequence = []
    for i in range(len(open_reading_frames)):
        if len(open_reading_frames[i]) > threshold:
            aminos = coding_strand_to_AA(open_reading_frames[i]) #return the amino_acids
            aa_sequence.append(aminos) #add it to the list
        i += 1
    return aa_sequence ## return the list containing the amino acid sequence encoded longer than treshold

if __name__ == "__main__":
    import doctest
    # Importing dna
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna)) #execute function
    # doctest.testmod(verbose=True)
    #doctest.run_docstring_examples(find_all_ORFs, globals(), verbose=True)
