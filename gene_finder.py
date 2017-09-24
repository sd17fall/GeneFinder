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
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'

    """
    complement_strand = ""
    for pair in nucleotide:
        if pair == 'A':
            complement_strand = complement_strand + 'T'
        if pair == 'C':
            complement_strand = complement_strand + 'G'
        if pair == 'T':
            complement_strand = complement_strand + 'A'
        if pair == 'G':
            complement_strand = complement_strand + 'C'
    #print (complement_strand)
    return(complement_strand)

    # TODO: implement this
    pass


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
    rev_complement =""
    complement = get_complement(dna)
    strand_length = len(complement)
    while (strand_length > 0):  #starts at the back of the strand
        rev_complement = rev_complement +  complement[strand_length-1]
        strand_length -= 1      #moves backwards
    return (rev_complement)
    # TODO: implement this
    pass


def rest_of_ORF(dna): # need to get this checked
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

    >>> rest_of_ORF("TAG")
    ''

    """
    #stop codons are TAA, TAG, TGA
    # TODO: implement this
    #stop_codons = ['TAA','TAG','TGA']
    stop_codons = ['TAA','TAG','TGA']

    for i in range (0, len(dna)-1,3): # move though the strand every three bacepairs

        if ( len(dna) - i >=3):       # only check for stop if there are three or more bace pairs left
            if dna[i:i+3] in stop_codons:   #stop if you find a stop codon

                return dna[0:i]         #return dna before the codon


    return dna
    pass


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
    >>> find_all_ORFs_oneframe("ATGATGCATGAATGTAGATAGATGTGCCC")
    ['ATGATGCATGAATGTAGA', 'ATGTGCCC']
    """
    all_pairs = len(dna)  #records the length of the total DNA strand
    total_checked = 0     #the number of strands checked.
    orfs = []             #list of found orfs
   
    i =0 #The number of pairs that have been checked
    
   
    while i <= all_pairs:  #while loop goes through the length of the strand
        stran = ""
        if dna[i:i+3] == "ATG": #if a start codon is found
            stran = "ATG"       # records the start codon
            #print (dna[i:])
            i = i + 3           #skipps the start codon
            stran = stran +  rest_of_ORF(dna[i:]) #add the rest of the orf to the start codon
            orfs.append(stran)  #saves the found orf
            i = i + 3   #skips stop codon
        i = i + 3       #moves to next framm
    return orfs         #resturns found orfs

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
    orf = []
    stran1 = find_all_ORFs_oneframe(dna)    #finds orfs in frame 1
    orf = (stran1)
    orf = orf + (find_all_ORFs_oneframe(dna[1:])) #finds orfs in frame 2
    orf = orf +(find_all_ORFs_oneframe(dna[2:])) #finds orfs in frame 3
    return orf          #returns found codon
    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    orf = []
    orf = orf +  find_all_ORFs(dna)  #finds orfs all refrence frames
    orf = orf + find_all_ORFs(get_reverse_complement(dna)) #finds orf on the other side of the stran
    return orf
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    return max(find_all_ORFs_both_strands(dna), key=len); #returns longest orf
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles

        """
    max_orf = ""
    temp = ""
    for i in range(0, num_trials):
        random_dna = shuffle_string(dna) #suffles DNA
        temp = longest_ORF(random_dna)   # finds the longest orf in the shuffled dna
        if len(temp) > len(max_orf):     #if the current orf is longer than the previose longest 
            max_orf = str(temp)          #make the current orf the longest



    return len(max_orf) #return the longest orf of all the orfs
    # TODO: implement this
    pass


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
    # >>> coding_strand_to_AA("ATGCGAATGTAGCATCAAA")
    # 'MPA'

    ans = ''
    # for codon in codons:
    #     for c in codon:
    #         for i in range(0,len(dna),3):
    #
    #             if dna[i:i+3] == c:
    #                 ans = ans + aa[codons.index(codon)]
    #
    #
    # return ans
    for i in range(0,len(dna),3): #amino acids are three base pairs long
        if len(dna[i:i+3]) < 3:   #read as long as there are atleast 3 bace pairs left
            continue
        ans = ans + aa_table[dna[i:i+3]]    #reads three bace pairs and matches it to an amino acid
    return                        #return amino acid
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    found_dna = []
    long_acids = []
    threshold = longest_ORF_noncoding(dna, 1500) # longest random orf after 1500 samples
    # print(threshold)

    #print(dna)
    #print(longest_ORF_noncoding(dna,1000)) the average length is 332 110 amino acids
    found_dna = find_all_ORFs_both_strands(dna) #finds all of the orfs in the dna
    # print(found_dna)
    found_acids = []
    for c in found_dna:
        found_acids.append(coding_strand_to_AA(c))  #translates bacepairs into amino acids
    for c in found_acids:
        if len(c) >= threshold/3:   #threshold was measured in bace pairs there are 3 base pairs per amino acid
            long_acids.append(c)    #list of amino acid orfs that are longer than the threshold
    print("threshold: " + str(threshold))   
    print ("Found Acids: ")
    
    return long_acids



if __name__ == "__main__":
    import doctest
    # print(len(dna))

    gene_finder(dna)
