################################################################################
###############################     WARNING!     ###############################
#### BOTH PYTHON FILE AND FASTA TXT FILE MUST BE INCLUDED IN SAME DIRECTORY ####
######  IF YOU COULD.... MAYBE IT'LL BE WORK WHEN IT'S NOT SAME DIRECTORY ######
################################################################################

import os
import sys
import time
from PrivateFunction import GlobalAlignment, Chunks

start = time.time()

# Exit code 2 : If the input file in the command-line argument does not exist
input_file = sys.argv[1]
if os.path.isfile(input_file) == False:
    sys.stderr.write('No input file.')
    exit(2)

# Read file after read availability is checked
f = open(input_file, 'r', encoding='cp949')
data = f.read()
f.close()

# Exit code 3 : If the input file has nothing
if len(data) == 0:
    sys.stderr.write('No protein sequence.')
    exit(3)

data_split = data.split('\n')
len_data_split = len(data_split)

FASTA_1 = []
FASTA_2 = []
FASTA_cnt = 0

# Start of seq2list loop
for i in range(0, len_data_split):
    if data_split[i][0] == '>':
        FASTA_cnt += 1
        # Break seq2list loop when read more than 2 sequences
        if FASTA_cnt > 2:
            break
        continue
    else:
        if FASTA_cnt == 1:
            FASTA_1 += data_split[i].upper()
        if FASTA_cnt == 2:
            FASTA_2 += data_split[i].upper()
# End of seq2list loop

# Exit code 4 : If input file does not follow the FASTA format
if FASTA_cnt == 0:
    sys.stderr.write('No correct format.')
    exit(4)

# Exit code 5 : If the input file has only one sequence
if FASTA_cnt == 1:
    sys.stderr.write('Need one more sequence.')
    exit(5)

# Do global pairwise sequence algorithm
score_mat = GlobalAlignment(FASTA_1, FASTA_2)

FASTA_1_len = len(FASTA_1) - 1
FASTA_2_len = len(FASTA_2) - 1
match = 1
mismatch = -1
gap = -1

rev_FASTA_1 = []
rev_FASTA_2 = []

# Find the REAL sequence of FASTA codes
while (FASTA_1_len >= 0 or FASTA_2_len >= 0):
    if (FASTA_1_len >= 0 and FASTA_2_len >= 0 and score_mat[FASTA_1_len][FASTA_2_len] == score_mat[FASTA_1_len - 1][FASTA_2_len - 1] + match):
        rev_FASTA_1.append(FASTA_1[FASTA_1_len])
        rev_FASTA_2.append(FASTA_2[FASTA_2_len])
        FASTA_1_len -= 1
        FASTA_2_len -= 1
    elif (FASTA_1_len >= 0 and FASTA_2_len >= 0 and score_mat[FASTA_1_len][FASTA_2_len] == score_mat[FASTA_1_len - 1][FASTA_2_len - 1] + mismatch):
        rev_FASTA_1.append(FASTA_1[FASTA_1_len])
        rev_FASTA_2.append(FASTA_2[FASTA_2_len])
        FASTA_1_len -= 1
        FASTA_2_len -= 1
    elif (FASTA_1_len >= 0 and score_mat[FASTA_1_len][FASTA_2_len] == score_mat[FASTA_1_len - 1][FASTA_2_len] + gap):
        rev_FASTA_1.append(FASTA_1[FASTA_1_len])
        rev_FASTA_2.append('-')
        FASTA_1_len -= 1
    else:
        rev_FASTA_1.append('-')
        rev_FASTA_2.append(FASTA_2[FASTA_2_len])
        FASTA_2_len -= 1

# Reverse since I saved them from end of list
seq_FASTA_1 = list(reversed(rev_FASTA_1))
seq_FASTA_2 = list(reversed(rev_FASTA_2))

# Print sequence per 60 letters
Chunks(seq_FASTA_1, seq_FASTA_2, 60)

# Print final score of global pairwise sequence alignment
print('Final score : ', score_mat[-1][-1])

print('\nTotal time : ', round(time.time() - start, 3), 'sec')
