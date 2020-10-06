def GlobalAlignment(FASTA_1, FASTA_2):
    match = 1
    mismatch = -1
    gap = -1

    FASTA_1_len = len(FASTA_1)
    FASTA_2_len = len(FASTA_2)

    # Make 2D list filled with 0
    score_mat = [[0 for i in range(FASTA_2_len + 1)] for j in range(FASTA_1_len + 1)]

    # Set first row and column with -1 added
    for i in range(1, FASTA_1_len + 1):
        score_mat[i][0] = score_mat[i-1][0] - 1
    for j in range(1, FASTA_2_len + 1):
        score_mat[0][j] = score_mat[0][j-1] - 1

    # Do the global pairwise sequence alignment algorithm
    for i in range(1, FASTA_1_len + 1):
        for j in range(1, FASTA_2_len + 1):
            score_row = score_mat[i][j-1] + gap
            score_col = score_mat[i-1][j] + gap

            if FASTA_1[i-1] == FASTA_2[j-1]:
                score_diag = score_mat[i-1][j-1] + match
                flag = 1
            else:
                score_diag = score_mat[i-1][j-1] + mismatch
                flag = 0

            score_list = [score_row, score_col, score_diag]
            score_mat[i][j] = max(score_list)

    return score_mat


def Chunks(list1, list2, n):
    for i in range(0, len(list1), n):
        print(''.join(list1[i:i+n]))
        print(''.join(list2[i:i+n]))
        print('\n')
