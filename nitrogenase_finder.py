"""
Extension on first project for Olin Software Design Fall 2017

@author: Emma Westerhoff

"""

from load import load_nitrogenase_seq, load_metagenome

def longest_common_substring(string1, string2): #s length r, t ength n
    """ Computes the longest common substring using dynamic programming

    >>> longest_common_substring('abcdefgqwertyuiop', 'xyabcdjipqwertyuiop')
    'gqwertyuiop'
    """

    x = len(string1)
    y = len(string2)

    L = [[None]*(y) for a in range(x)]

    z = 0
    ret = ''

    for i in range(0, x):
        for j in range(0, y):
            if string1[i] == string2[j]:
                if i == 0 or j == 0:
                    L[i][j] = 0
                else:
                    L[i][j] = L[i-1][j-1] + 1
                if L[i][j] > z:
                    z = L[i][j]
                    ret = string1[i-z:i+1]
                #elif L[i][j] == z:
                    #ret.append(string1[i-z:i+1])
            else:
                L[i][j] = 0
    return ret

def nitrogen_fixation(x):
    #TODO: implement this
    pass

if __name__ == "__main__":
    nitrogenase = load_nitrogenase_seq()
    metagenome = load_metagenome()
    #metagenome is of form [('some info', 'actual sequence')]
    #transform metagenome to proper form
    
    #import doctest
    #doctest.run_docstring_examples(longest_common_substring, globals(), verbose=True)
    #longest = longest_common_substring(nitrogenase, metagenome)
    #print(longest)
