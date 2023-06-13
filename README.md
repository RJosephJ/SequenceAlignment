# SequenceAlignment
Sequence Alignment Code using Needleman-Wunsch Algorithm
"""
    @Ricardo Rivera Sanchez
    Needleman-Wunsch Algorithm for Global Alignment
    
    This is an implementation of the Needleman-Wunsch algorith for global alignment. This program
    takes a csv file as input from the command line containing two sequences that will be used for the alignment. The program 
    works by creating two matrices to be filled for the alignment. One matrix is the standard main_matrix which will be filled
    out with the Needleman-Wunsch algorithm which uses a system of matches and mismatches. For this implementation our mismatch will       have a penalty of -1,
    a match a reward of 1 and a gap penalty of -2.The second matrix is filled with the matches and mismatches established before. For     the backtracking section of
    the alignment we check the values in reverse, tracing back to the beginning. For this the program has established a left bias 
    so that in the case of the same value coming from the diagonal, left or upwards we prioritize the left first, then the upwards and     lastly the diagonal. If the value at the given iteration in the matrix that we are checking has the same value as the sum plus a       mismatch penalty to the left, then the traceback will move to the left. If
    the value of the iteration is the same as the sum of the value of the upward position of the matrix and its mismatch value, the       traceback will move
    upwards. Lastly, if the iteration we are on in the matrix is the sum of the diagonal and the reward value/penalty then the        traceback will move diagonally. 
    There were cases where all of the letters of the sequence where not checked, because the matrix at the beginning either moved to the right or downwards when finding its maximum instead of the
    usual diagonal movement in the beginning. For these special cases there are two iterations created depending on the case. 
    If the the first sequence had unchecked letters at the end of the traceback, we add a dash to the second alignment for how many unchecked letters there were. 
    The same applies but in reverse when the second sequence had unchecked letters due to the behaviour of the matrix.
    

"""
