"""
This module maps codons to amino acids.

Attributes:
    aa_table (dict(str, str)):
        A map from codons to amino acids. Codons are represented as strings
        in IUPAC notation, e.g. 'TTT'. Amino acids are represented by their
        single-letter IUPAC nucleobase abbreviations, e.g. 'F' for
        phenylalanine.

        Examples:

        >>> aa_table['TTC']
        'F'
        >>> aa_table['TTA']
        'L'
        >>> aa_table['GGA']
        'G'

References:
    * DNA codon table: https://en.wikipedia.org/wiki/DNA_codon_table
    * IUPAC notation: https://en.wikipedia.org/wiki/Nucleic_acid_notation

"""

# amino acids, in the same order as in `codons`. This is used to construct
# `aa_table`.
amino_acids = [
    'F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y', '|', 'H', 'Q', 'N', 'K',
    'D', 'E', 'C', 'W', 'R', 'G']

# A list of lists of codons, in the same order as `amino_acids`.
codons = [
    ['TTT', 'TTC'],
    ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    ['ATT', 'ATC', 'ATA'],
    ['ATG'],
    ['GTT', 'GTC', 'GTA', 'GTG'],
    ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    ['CCT', 'CCC', 'CCA', 'CCG'],
    ['ACT', 'ACC', 'ACA', 'ACG'],
    ['GCT', 'GCC', 'GCA', 'GCG'],
    ['TAT', 'TAC'],
    ['TAA', 'TAG', 'TGA'],
    ['CAT', 'CAC'],
    ['CAA', 'CAG'],
    ['AAT', 'AAC'],
    ['AAA', 'AAG'],
    ['GAT', 'GAC'],
    ['GAA', 'GAG'],
    ['TGT', 'TGC'],
    ['TGG'],
    ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    ['GGT', 'GGC', 'GGA', 'GGG']
]

# create a dictionary lookup table for mapping codons into amino acids
# See Think Python Ch. 11
# http://greenteapress.com/thinkpython/html/thinkpython012.html
# http://greenteapress.com/thinkpython/html/thinkpython012.html
aa_table = {}
for i in range(len(amino_acids)):
    for codon in codons[i]:
        aa_table[codon] = amino_acids[i]

# Use `enumerate` to avoid having to compute the size of the list,
# and then select each of its items as a separate operation.
# `enumerate` produces a sequence of (index, item) pairs that
# `i` and `nucleotide_codons` bind to.
# See: https://docs.python.org/3/library/functions.html#enumerate
aa_table = {}
for i, nucleotide_codons in enumerate(codons):
    for codon in nucleotide_codons:
        aa_table[codon] = amino_acids[i]

# `zip` produces a sequence of pairs: a pair that contains the first item from
# each of `codons` and `amino_acids`, then the second item from each, and so
# on.
# See: https://docs.python.org/3/library/functions.html#zip
aa_table = {}
for nucleotide_codons, amino_acid in zip(codons, amino_acids):
    for codon in nucleotide_codons:
        aa_table[codon] = amino_acid

# instead of creating an empty dictionary and then filling it, use a
# dictionary comprehension to create the dictionary in one swell foop.
aa_table = {codon: aa
            for aa_codons, aa in zip(codons, amino_acids)
            for codon in aa_codons}

# instead of including `codons` and `amino_acids` as above, read
# the map from a CSV ("commma-separated values") file.
# The file looks like this:
#   AAA,K
#   AAC,N
#   AAG,K
# etc.
aa_table = dict(line.strip().split(',')
                for line in open('codon_nucleotides.csv').readlines())



def fn(a, b):
    if test():
        return


#!/usr/bin/env python3


# The Pandas library can read CSV and Excel tables, and manipulate tables
# within a Python program. It's overkill for this purpose, but
import doctest

# very powerful in general.
# See: http://pandas.pydata.org
import pandas as pd

aa_table = pd.read_csv('codon_nucleotides.csv', header=None, index_col=0)[1]

# The previous solution, unlike all the others, doesn't actually create
# a Python dict. It creates a (Pandas) Series, which can be used the same
aa_table = dict(pd.read_csv('codon_nucleotides.csv',
                            header=None, index_col=0)[1])
# do this instead:
aa_table = dict(pd.read_csv('codon_nucleotides.csv',
                            header=None, index_col=0)[1])

# Finally, we could just list out the key-value pairs in the dictionary.
aa_table = {
    'AAA': 'K',
    'AAC': 'N',
    'AAG': 'K',
    'AAT': 'N',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'AGA': 'R',
    'AGC': 'S',
    'AGG': 'R',
    'AGT': 'S',
    'ATA': 'I',
    'ATC': 'I',
    'ATG': 'M',
    'ATT': 'I',
    'CAA': 'Q',
    'CAC': 'H',
    'CAG': 'Q',
    'CAT': 'H',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'GAA': 'E',
    'GAC': 'D',
    'GAG': 'E',
    'GAT': 'D',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'TAA': '|',
    'TAC': 'Y',
    'TAG': '|',
    'TAT': 'Y',
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TGA': '|',
    'TGC': 'C',
    'TGG': 'W',
    'TGT': 'C',
    'TTA': 'L',
}
    'TTG': 'L',
    'TTT': 'F'
}


doctest.testmod()
