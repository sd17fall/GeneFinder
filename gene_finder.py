# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Aditya Kaushika

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

def count_v1(dna, base):
    dna = list(dna)  # convert string to list of letters
    i = 0            # counter
    for c in dna:
        if c == base:
            i += 1
    return i


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s))) # allows us to shuffle the letters in the code

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

def get_complement(nucleotide):
    """ Returns the complementary nucleotide nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    for letter in nucleotide: #Go through every letter
        if letter == 'A': #If Letter is A:
            return 'T'      #Then complement it with T
        elif letter == 'T': #Repeat with different letters.
            return 'A'
        elif letter == 'G':
            return 'C'
        else:
            return 'G'

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
    Dna_reversed = dna[::-1]   #Start counting backwards
    d=''   #Create new sting
    for letter in Dna_reversed:  #Look at each backwards letter
        reverse = get_complement(letter) #Get the complement of that letter
        d+= reverse   #Add it to our new list
    return d   #Visualize our list


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
    for i in range(0,len(dna),3): #look at the range from the 0th place to the 2nd place
        if (dna[i:i+3]== "TAA") or (dna[i:i+3]== "TAG") or (dna[i:i+3]== "TGA") : #Identify the stop codons
            return dna[0:i] #Stop if you see a stop codon
    return dna #Display results


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
    a=[]  #Create List
    rol=0  #start our remaining letters at 0
    for i in range(0,len(dna),3): # Look at three letters at a time
        codon= dna[i:i+3] #assign value to codon
        if (codon =='ATG' and rol<=0): #Restrictions for our if statement
            orf= rest_of_ORF(dna[i:]) #assign value to orf
            a.append(orf)             #add on to our list
            rol=len(orf)              #add letters to our remaining letters
        if  (rol > 0):                #ask if there are any remaining letters
            rol = rol - 3             #if there are, take three away
    return a                           #Display results



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
    a=[]     #New List
    for i in range(3): #Look at all three frams, instead of one
        a.extend(find_all_ORFs_oneframe(dna[i:])) #Add onto list a

    # for i in range(0,len(dna),3):
    #     if (dna[i:i+3]=='ATG'):
    #         a.append(rest_of_ORF(dna[i:]))
    #     elif (dna[i+1:i+4]=='ATG'):
    #         a.append(rest_of_ORF(dna[i:]))
    #     elif (dna[i+2:i+5]=='ATG'):
    #         a.append(rest_of_ORF(dna[i:]))
    return a   #Display results


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    a=[]                                        #New List
    a.extend(find_all_ORFs(dna))                #add last functions results to new list
    dna_reversed = dna[::-1]                    #Reverse direction of letters
    d=""                                        #New String
    for letter in dna_reversed:                 #look at each letter in new variable
        reverse = get_complement(letter)        #get the complement of reversed letters
        d+= reverse                             #add onto the string d
    a.extend(find_all_ORFs(d))                  #add onto the list a
    return a                                    #display results

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    count=0                                     #Start Count at 0
    for i in find_all_ORFs_both_strands(dna):   #For items in the list created above
        if len(i)>count:                        #Check their length to the last counted string
            count = len (i)                     #If it is the longest, change it
            a = i                               #make a the largest string
    return a                                    #Display a

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        """
    List_of_Longest_Orfs=[] #New List
    count = 0 #Start the count at 0
    for i in range(num_trials): #Tells the for statement how many times to run
        shuffled_dna = longest_ORF(shuffle_string(dna[i:])) #shuffle the DNA string
        List_of_Longest_Orfs.append(shuffled_dna)#add the shuffled DNA into a list
    for i in (List_of_Longest_Orfs): #function to look into the afformentioned list
        if len(i)>count: #Compare the length to the count of the previous word
            count = len (i) #make the longer length = to the new count
            List_of_Longest_Orfs = i #Tell the function where to get its list
    return count #Display the results for us to see


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
    a="" #New List for results to be put into
    for i in range(0,len(dna),3): #Tell the function how long to look for
        codon = dna[i:i+3] #Tell the function where to look
        amino_acid = aa_table[codon] #assign value to amino_acid
        a = a + amino_acid #continue the string
    return a #actually return the amino acid

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    a=[]  #creating an empty list
    threshold = longest_ORF_noncoding(dna,1500) #Assign Value to threshold
    Long_Orfs = len(longest_ORF(dna)) #assign value to Long_Orfs
    if Long_Orfs>threshold): #Compare values
        a.append(coding_strand_to_AA(dna)) #add to the list
    dna = load_seq("./data/X73525.fa") #obtaining genes
    print gene_finder(dna) #showing the list of Amino Acids


if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA, globals(),verbose=True)
    # doctest.
