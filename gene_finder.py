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
    """ Checks for each letter and returns the complement
        (might not be the most effiicent way) and returns
        'This is not a nucleotide' if any other letters are implented"""
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    pass


def get_reverse_complement(dna):
    reverse_dna = ''
     """ base string to add complements"""
    index = 0
     """ counter for while loop"""
    while index < len(dna):
        reverse_dna = reverse_dna + get_complement(dna[index])
        index = index + 1
    """ while loop that processes each letter of the string and applies the
        previous function get_complement() to add on to the newly made string
        reverse_dna"""
    return reverse_dna



    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    pass

#This one code took me an hour to write. jesus.
def rest_of_ORF(dna):
  index = 0
  newstr = ''
  stop_codon = ['TAG','TAA','TGA']
  """index is set as the counter, newstr is what the output is set to be, and
     the stop_codon variable is just for me to check for stop codons"""
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

    """The rest of the function checks every three codons and tries to match them
    to the stop codon. If it isn't a stop codon then it will add the letter to the string"""

        """ Takes a DNA sequence that is assumed to begin with a start
            codon and returns the sequence up to but not including the
            first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        TAG TAA and TGA


        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
    pass


""" V1
def find_all_ORFs_oneframe(dna):
  newstr1 = []
  newstr1.append(rest_of_ORF(dna))
  stop_codon = ['TAG','TAA','TGA']
  index = 0
  #newstr1 adds the first chunk to the list
  while dna[index:index+3] not in stop_codon:
    index = index + 3
  newstr1.append(dna[index+3:])
  #checks to find
  return newstr1
"""

#V2 all_ORFS_oneframe
def find_all_ORFs_oneframe(dna):
  newstr1 = []
  newstr1.append(rest_of_ORF(dna))
  stop_codon = ['TAG','TAA','TGA']
  index = 0
  #newstr1 adds the first chunk to the list
  while dna[index:index+3] not in stop_codon and (index < len(dna)-3):
      index = index + 3
  newstr1.append(dna[index+3:])
  return newstr1

#print(find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC"))
#print(find_all_ORFs_oneframe('AAAAAAAATAGAAA'))
#print(find_all_ORFs_oneframe('AAAATTTTTTAATTTTTTATTATA'))

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
    # TODO: implement this
    pass


#this code doesn't work yet. Still going through testing.
def cutlist(dna):
  list = dna
  index = 1
  while 'ATG' != dna[index:index + 3]:
    index = index + 1
  list = dna[index:]
  return list

print (cutlist("ATGCATGAATGTAG"))


def find_all_ORFs(dna):
  finallist = []
  list = dna
  while 'ATG' in list:
    if 'ATG' == list[0:3]:
      finallist.append(list)
      list = cutlist(list)
    else:
      list = cutlist(list)
  return finallist


print(find_all_ORFs("ATGCATGAATGTAG"))

print(find_all_ORFs("ATGCATGAATGTAG"))

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
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
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
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
