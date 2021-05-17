''' Assignment 2 - solving the Shortest Common Super-string (SCS) problem'''
'''by Pietro Zafferani'''

import random
import itertools

'''Creates randomly a DNA string of the desired length.'''


def getGenome(length: int) -> str:
    genome = ''.join(random.choice('ACTG') for i in range(length))
    return genome


'''Given a DNA stretch (genome) this functions creates all its substrings of a given length.'''


def getSubstrings(genome: str, length: int) -> list:
    l = []
    # define the number of substrings
    for i in range(len(genome) - length + 1):
        # add each substring to the final list
        l.append(genome[i: i + length])
    return l


'''This function exploit the itertools module to generate all the permutations of a given list of sub-strings'''


def getPermutations(Set: list) -> list:
    for perm in itertools.permutations(Set):
        # generators return a single permutation each function calling, so the PC memory usage is lighter
        yield list(perm)


'''This function checks whether a pair of sub-strings has a suffix/prefix type of overlap.'''


# seq1 -> check suffix
# seq2 -> check prefix
def overlap(seq1: str, seq2: str) -> int:
    best = 0
    # execute until the index falls inside the shortest substrings
    for i in range(1, min(len(seq1), len(seq2)) + 1):
        # check the presence of suffix/prefix overlap
        if seq1[-i:] == seq2[:i]:
            # this index indicates the length of the largest overlap
            best = i
    # return only the largest suffix/prefix overlap found
    return best


'''Creates the skeleton of the overlap matrix, made of only the first raw and first column. 
    Logically the matrix rows are equal to the columns.'''


def emptyMatrix(L: list) -> list:
    M = [['/']]
    for column in L:
        # creates first row's labels
        M[0].append(column)
        # creates the first element of each row
        row = [column]
        # append each row to the matrix
        M.append(row)
    return M


'''This function fills the overlap matrix, every cell belonging to 2 different sub-strings contains the overlapping
    score of the 2 sub-strings. The output is formed by several pieces of information that will be exploited
    by the next functions in order to accomplish the task.'''


def OverlapMatrix(substrings: list) -> list:
    
    # tells the merged string where to be inserted in the new set of strings
    index_Y = 0
    # defines the 2 strings to be merged in the Greedy function
    best_match = ('', '')
    # defines the best overlap-score
    best_score = 0
    # create the overlap matrix given a set of sub-strings
    OverlapM = emptyMatrix(substrings)
    # define the first row of the matrix
    first_row = OverlapM[0]
    
    # iterate over each row of the matrix
    for r in range(1, len(OverlapM)):
        row = OverlapM[r]
        # iterate over each column of each row
        for c in range(1, len(first_row)):
            column = OverlapM[0][c]
            # the diagonal of the overlap matrix must not be computed
            if r == c:
                # the self-overlap of a sub-string is set negative
                score=-1
                row.append(score)
            else:
                # compute the overlap score that is inserted in the matrix's cell
                score = overlap(row[0], column)
                row.append(score)
                # the best overlap score is stored, together with the respective pair of substring
                # if there are more than one highest-overlap score, the first one encountered is returned
            if score > best_score:
                best_score = score
                best_match = (row[0], column)
                # store the position in the set for the first string in the pair (converted becomes shifted by one)
                index_Y = r-1
                
    # return all the useful data retrieved
    return [OverlapM, best_score, best_match, index_Y]


'''It takes 2 substrings and their overlap score as arguments, then it merges them according to such score.'''


# seq1 provides the suffix
# seq2 provides the prefix
def merge(seq1: str, seq2: str, score) -> str:
    # compute the first untouched part of the upperstring
    firstpart = len(seq1) - score

    # return the merged upperstring
    return seq1[:firstpart] + seq2


'''This function is useful for testing the correct structure of the overlap matrix.'''


# For instance type: showMatrix(OverlapMatrix(Set_of_substrings))
def showMatrix(M: list) -> print:
    for i in M[0]:
        print(i)


'''This function applies a greedy strategy for the assembly of the SCS, however since this approach is heuristic,
    it is not sure that the best solution will always be returned. It implements recursion to produce the last 
    string remained.'''


def GreedySCS(Set: list) -> str:
    # base case of recursion
    if len(Set) == 1:
        return Set[0]

    # retrieve the output given by OverlapMatrix()
    RecursiveMatrix = OverlapMatrix(Set)
    seq1 = RecursiveMatrix[2][0]
    seq2 = RecursiveMatrix[2][1]
    score = RecursiveMatrix[1]
    oldSeq1 = RecursiveMatrix[3]

    # create the upperstring formed by the 2 sub-strings with highest overlap score
    upperString = merge(seq1, seq2, score)

    

    # remove the 2 sub-strings not present anymore
    Set.remove(seq1)
    Set.remove(seq2)
    
    # insert the merged string in the original set of sub-strings
    Set.insert(oldSeq1, upperString)
   
    # execute recursion until the base case is satisfied
    return GreedySCS(Set)


'''This function takes a list of sub-strings as argument and sequentially merges/concatenates them following the
    original order. Eventually a single upperstring is returned.'''


def serialMerger(L: list) -> str:
    # base case of recursion
    if len(L) == 1:
        return L[0]

    # merge the first 2 sub-strings of the set
    upperString = merge(L[0], L[1], overlap(L[0], L[1]))
    # remove the previous pair of sequences from the new set and add the merged one in their place
    new_L = L[2:]
    new_L.insert(0, upperString)

    # recursive case executed until only one string remains in the set
    return serialMerger(new_L)


'''This function applies a brute-force strategy to return the best SCS possible. It goes through all the possible 
    permutations of the original set of sub-strings, it creates a SCS for each of them and then returns only the
    shortest one.'''


def BruteForceSCS(Set: list) -> str:
    # the longest upperstring possible is given by the concatenation of all the substrings
    # set by default
    best_len = len(Set) * len(Set[0])
    # set by default as empty
    best_scs = ''
    # iterate over all the possible permutations
    for perm in getPermutations(Set):
        # compute the ordered merging/concatenating of each set
        outcome = serialMerger(perm)

        # check whether the new merged string is the shortest so far
        if len(outcome) < best_len:
            # store the new SCS and its length
            best_len = len(outcome)
            best_scs = outcome

    # return SCS after checking all possible upperstrings
    return best_scs


'''This function is the interface with the user, it asks if the user wants to test the algorithms a randomly
    generated set of DNA sub-strings or to insert a personal one.
    Then it asks which strategy of assembly the user wants to apply: Greedy or Brute.
    Eventually the SCS is returned.'''


def MAIN():
    question = input('Do you want to use a randomly generated set of substrings?  [yes, no] ')
    if str(question.lower()) == 'yes':
        question2 = input('Please insert genome\'s length: ')
        strings = getSubstrings(getGenome(int(question2)), int(question2)//3)
        print(strings)
    else:
        strings = eval(input('''Please insert your set of substrings in this format: [\'x\',...,\'z\'] '''))
    approach = input('Select assembling approach [Greedy, Brute]: ')
    print()
    if str(approach) == 'Greedy':
        print(GreedySCS(strings))
    else:
        print('The process  may require some time...')
        print(BruteForceSCS(strings))


if __name__ == '__main__':
    #show a random genome of lenght n -> print(getGenome(n))

    #show set of k-mers -> print(getSubstrings(genome, k))

    #print the overlap matrix -> showMatrix(OverlapMAtrix(substrings))

    #assembly the SCS with a greedy strategy -> print(GreedySCS(substrings))

    #assembly the SCS with a Brute force strategy -> print(BruteForceSCS(substrings))

    # instances of sub-strings sets
    # ['AAA','AAB','ABB','BBB','BBA']
    # ['BBB', 'AAA', 'ABB']
    # ['AAA','ABB','BBB','AAB','BBA']
    # ['CGATTT', 'GATTTT', 'ATTTTC', 'TTTTCT', 'TTTCTT']
   
    # run the file
    MAIN()
