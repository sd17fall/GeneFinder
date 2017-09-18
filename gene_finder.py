"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffle the characters in the input string.

    NOTE: this is a helper function, you do not
    have to modify this in any way."""
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """Return the complementary nucleotide.

    Args:
        nucleotide: a nucleotide (A, C, G, or T) represented as a string

    Returns:
        The complementary nucleotide

    Examples:
        >>> get_complement('A')
        'T'
        >>> get_complement('C')
        'G'
    """
    # TODO: implement this
    pass


def get_reverse_complement(dna):
    """Compute the reverse complementary sequence of DNA for the specified DNA sequence.

    Args:
        dna: a DNA sequence represented as a string

    Returns:
        The reverse complementary DNA sequence represented as a string.

    Examples:
        >>> get_reverse_complement("ATGCCCGCTTT")
        'AAAGCGGGCAT'
        >>> get_reverse_complement("CCGCGTTCA")
        'TGAACGCGG'
    """
    # TODO: implement this
    pass


def rest_of_ORF(dna):
    """Given a DNA sequence that begins with a start codon, return the ORF.

    Takes a DNA sequence that is assumed to begin with a start
    codon and returns the sequence up to but not including the
    first in frame stop codon.  If there is no in frame stop codon,
    returns the whole string.

    Args:
        dna: a DNA sequence

    Returns:
        The open reading frame represented as a string

    Examples:
        >>> rest_of_ORF("ATGTGAA")
        'ATG'
        >>> rest_of_ORF("ATGAGATAGG")
        'ATGAGA'
    """
    # TODO: implement this
    pass


def find_all_ORFs_oneframe(dna):
    """Find non-nested ORFs in the default frame of the sequence.

    Finds all non-nested open reading frames in the given DNA
    sequence and returns them as a list.  This function should
    only find ORFs that are in the default frame of the sequence
    (i.e. they start on indices that are multiples of 3).
    By non-nested we mean that if an ORF occurs entirely within
    another ORF, it should not be included in the returned list of ORFs.

    Args:
        dna: a DNA sequence

    Returns:
        A list of non-nested ORFs

    Examples:
        >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
        ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this
    pass


def find_all_ORFs(dna):
    """Find non-nested ORFs in all frames of the sequence.

    Finds all non-nested open reading frames in the given DNA sequence in
    all 3 possible frames and returns them as a list.  By non-nested we
    mean that if an ORF occurs entirely within another ORF and they are
    both in the same frame, it should not be included in the returned list
    of ORFs.

    Args:
        dna: a DNA sequence

    Returns:
        A list of non-nested ORFs

    Examples:
        >>> find_all_ORFs("ATGCATGAATGTAG")
        ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """Find all non-nested ORFs in both strands.

    Find all non-nested open reading frames in the given DNA sequence on both strands.

    Args:
        dna: a DNA sequence

    Returns:
        A list of non-nested ORFs

    Examples:
        >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
        ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass


def longest_ORF(dna):
    """Find the longest ORF on both strands of the specified DNA.

    Returns:
        The longest ORF, as a string

    Examples:
        >>> longest_ORF("ATGCGAATGTAGCATCAAA")
        'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """Computer the maximum length of the longest ORF from a set of shuffled sequenced.

    Compute the maximum length of the longest ORF over num_trials shuffles
    of the specified DNA sequence.

    Args:
        dna: a DNA sequence
        num_trials: the number of random shuffles

    Returns:
        The maximum length longest ORF
    """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """Compute the protein encoded by a sequence of DNA.

    Computes the protein encoded by a sequence of DNA.  This function
    does not check for start and stop codons (it assumes that the input
    DNA sequence represents an protein coding region).

    Args:
        dna: a DNA sequence represented as a string

    Returns:
        A string containing the sequence of amino acids encoded by the
        the input DNA fragment.

    Examples:
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """Returns the amino acid sequences that are likely coded by the specified sequence.

    Args:
        dna: a DNA sequence

    Returns:
        A list of all amino acid sequences coded by the sequence DNA.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
