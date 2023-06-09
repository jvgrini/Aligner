import numpy as np
from Bio import SeqIO

#Smith-Waterman alignment using constant match/mismatch values

def local_align(seq_file, m_score = 2, mm_score = -1, g_score = -1):

    seq_1 = ''
    seq_2 = ''
    seq_1_desc = ''
    seq_2_desc = ''

    filename = seq_file
    numseqs = 0
    for record in SeqIO.parse(filename, 'fasta'):
        numseqs += 1

    count = 0
    if numseqs == 2:
        for record in SeqIO.parse(filename, 'fasta'):
            if count == 0:
                seq_1 = str(record.seq)
                seq_1_desc = str(record.id)
            else:
                seq_2 = str(record.seq)
                seq_2_desc = str(record.id)
            count+=1


    m1 = np.zeros((len(seq_1)+1, len(seq_2)+1))

    #Fill matrix
    max_score = 0
    max_pos = None
    for i in range(1,len(seq_1)+1):
        for j in range(1, len(seq_2)+1):
            if seq_1[i-1] == seq_2[j-1]:
                score = m_score
            else:
                score = mm_score

            m1[i][j] = max(
                0,
                m1[i-1][j-1] + score,
                m1[i-1][j] + g_score,
                m1[i][j-1] + g_score
            )

            if m1[i][j] > max_score:
                max_score = m1[i][j]
                max_pos = (i, j)         
    a1 = ''
    a2 = ''

    i, j = max_pos
    its1 = 0
    its2 = 0
    while m1[i][j] != 0:
        if m1[i][j] == m1[i-1][j] + g_score:
            a1 = seq_1[i-1] + a1
            a2 = "-" + a2
            i -= 1
            its1 +=1
        elif m1[i][j] == m1[i][j-1] + g_score:
            a1 = "-" + a1
            a2 = seq_2[j-1] + a2
            j -= 1
            its2 +=1
        else:
            a1 = seq_1[i-1] + a1
            a2 = seq_2[j-1] + a2
            i -= 1
            j -= 1
            its1 +=1
            its2 +=1
        
    start_seq1 = max_pos[0] - its1 +1
    start_seq2 = max_pos[1] - its2 +1
        #gen line
    matchlines = ''
    for j in range(len(a1)):
        if a1[j] == a2[j]:
            matchlines = matchlines + '|'
        else:
            matchlines = matchlines + ' '
        j +=1

     #Score alignment
    identity = 0
    positives = 0
    gaps = 0
    for i in range(len(matchlines)):
        if matchlines[i] == '|':
            identity +=1
        if matchlines[i] == '+':
            positives +=1
    for i in range(len(a1)):
        if a1[i]== '-':
            gaps +=1
    for i in range(len(a2)):
        if a2[i]== '-':
            gaps+=1
    
    positives = positives + identity

    return a1, a2, matchlines, seq_1_desc, seq_2_desc, start_seq1, start_seq2, identity, positives, gaps
