import csv
import numpy as np

"""
    @Ricardo Rivera Sanchez
    Needleman-Wunsch Algorithm for Global Alignment
    
    
    This is an implementation of the Needleman-Wunsch algorith for global alignment. This program
    takes a csv file as input from the command line containing two sequences that will be used for the alignment. The program 
    works by creating two matrices to be filled for the alignment. One matrix is the standard main_matrix which will be filled
    out with the Needleman-Wunsch algorithm which uses a system of matches and mismatches. For this implementation our mismatch will have a penalty of -1,
    a match a reward of 1 and a gap penalty of -2.The second matrix is filled with the matches and mismatches established before. For the backtracking section of
    the alignment we check the values in reverse, tracing back to the beginning. For this the program has established a left bias 
    so that in the case of the same value coming from the diagonal, left or upwards we prioritize the left first, then the upwards and lastly the
    diagonal. If the value at the given iteration in the matrix that we are checking has the same value as the sum plus a mismatch penalty to the left, then the traceback will move to the left. If
    the value of the iteration is the same as the sum of the value of the upward position of the matrix and its mismatch value, the traceback will move
    upwards. Lastly, if the iteration we are on in the matrix is the sum of the diagonal and the reward value/penalty then the traceback will move diagonally. 
    There were cases where all of the letters of the sequence where not checked, because the matrix at the beginning either moved to the right or downwards when finding its maximum instead of the
    usual diagonal movement in the beginning. For these special cases there are two iterations created depending on the case. 
    If the the first sequence had unchecked letters at the end of the traceback, we add a dash to the second alignment for how many unchecked letters there were. 
    The same applies but in reverse when the second sequence had unchecked letters due to the behaviour of the matrix.
    

"""

def needle(filename):

    #Open csv file for extracting the sequences to be used in the algorithm
    with open(filename, 'r') as csv_file:
        csv_file = csv.reader(csv_file)
        for line in csv_file:
            if line[0].isupper():
                #Store both sequences on two variables for ease of use
                sequence_1 = line[0]
                sequence_2 = line[1]

                main_matrix = np.zeros((len(sequence_1) + 1, len(sequence_2) + 1))
                match_checker_matrix = np.zeros((len(sequence_1), len(sequence_2)))
                match_reward = 1
                mismatch_penalty = -1
                gap_penalty = -2
                for i in range(len(sequence_1)):
                    for j in range(len(sequence_2)):
                        if sequence_1[i] == sequence_2[j]:
                            match_checker_matrix[i][j] = match_reward
                        else:
                            match_checker_matrix[i][j] = mismatch_penalty
                #Initialisation of matrix
                for i in range(len(sequence_1) + 1):
                    main_matrix[i][0] = i * gap_penalty
                for j in range(len(sequence_2) + 1):
                    main_matrix[0][j] = j * gap_penalty

                # Matrix Filling
                for i in range(1, len(sequence_1) + 1):
                    for j in range(1, len(sequence_2) + 1):
                        main_matrix[i][j] = max(main_matrix[i - 1][j - 1] + match_checker_matrix[i - 1][j - 1],main_matrix[i - 1][j] + gap_penalty,main_matrix[i][j - 1] + gap_penalty)


                #Traceback
                aligned_1 = ""
                aligned_2 = ""

                ti = len(sequence_1)
                tj = len(sequence_2)
                #For each length of the sequences check the matrix in reverse to catch the traceback

                #By checking if the current value of the matrix in the present iteration is the sum of the value of the
                #diagonal, left or upwards plus the penalty/reward, then we can safely move to that direction that gave us the number.
                #It works like an arrow pointing at the direction it came from
                while (ti > 0 and tj > 0):
                    if (tj > 0 and main_matrix[ti][tj] == main_matrix[ti][tj - 1] + gap_penalty):

                        aligned_1 = "-" + aligned_1
                        aligned_2 = sequence_2[tj - 1] + aligned_2

                        tj = tj - 1

                    elif (ti > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj] + gap_penalty):

                        aligned_1 = sequence_1[ti - 1] + aligned_1
                        aligned_2 = "-" + aligned_2

                        ti = ti - 1






                    elif (ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj - 1] +
                          match_checker_matrix[ti - 1][tj - 1]):

                        aligned_1 = sequence_1[ti - 1] + aligned_1
                        aligned_2 = sequence_2[tj - 1] + aligned_2

                        ti = ti - 1
                        tj = tj - 1
                #Special cases for leftover letters in the alignment
                #In these special cases where the matrix moves to the right or downwards at the beginning of the sequence
                #The len of the sequences are not 0
                #For each unchecked letter, a "-" is added depending on the sequence that was not checked completely
                if (ti > 0):
                    while (ti > 0):
                        aligned_2 = "-" + aligned_2
                        aligned_1 = sequence_1[ti - 1] + aligned_1
                        ti = ti - 1
                if (tj > 0):
                    while (tj > 0):
                        aligned_1 = "-" + aligned_1
                        aligned_2 = sequence_2[tj - 1] + aligned_2
                        tj = tj - 1

                print(aligned_1, aligned_2, int(main_matrix[len(sequence_1)][len(sequence_2)]))








