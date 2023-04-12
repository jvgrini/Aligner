import numpy as np
from Bio import SeqIO
import blosum as bl

#Needleman-Wunch alignment with constant match/mismatch values


def global_align(seq_file, m_score = 1, mm_score = -1, g_score = -2):
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
#create matrix
    m1 = np.zeros((len(seq_1)+1, len(seq_2)+1))
    m2 = np.zeros((len(seq_1), len(seq_2)))


    for i in range(len(seq_1)):
        for j in range(len(seq_2)):
            if seq_1[i] == seq_2[j]:
                m2[i][j] = m_score
            else:
                m2[i][j] = mm_score

    #Fill matrix

    for i in range(len(seq_1)+1):
        m1[i][0] = i*g_score
    for i in range(len(seq_2)+1):
        m1[0][i] = i*g_score

    for i in range(1,len(seq_1)+1):
        for j in range(1, len(seq_2)+1):
            m1[i][j] = max(m1[i-1][j-1]+m2[i-1][j-1],m1[i-1][j]+g_score, m1[i][j-1]+g_score)

    #find alignment
    a1 = ''
    a2 = ''

    l1 = len(seq_1)
    l2 = len(seq_2)

    while(l1 > 0 and l2 > 0):
        if(l1 > 0 and l2 > 0 and m1[l1][l2] == m1[l1-1][l2-1] + m2[l1-1][l2-1]):
            a1 = seq_1[l1-1] + a1
            a2 = seq_2[l2-1] + a2

            l1 = l1 -1
            l2 = l2 -1
        elif(l1>0 and m1[l1][l2] == m1[l1-1][l2] + g_score):
            a1 = seq_1[l1-1] + a1
            a2 = "-" + a2

            l1 = l1-1
        else:
            a1 = '-' + a1
            a2 = seq_2[l2-1] + a2

            l2 = l2-1

    #gen line
    matchlines = ''
    for i in range(len(a1)):
        if a1[i] == a2[i]:
            matchlines = matchlines + '|'
        else:
            matchlines = matchlines + ' '
        i +=1

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

    

    return a1, a2, matchlines, seq_1_desc, seq_2_desc, positives, identity, gaps