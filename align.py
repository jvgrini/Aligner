from Bio import SeqIO
import numpy as np
import argparse
from pathlib import Path

from scripts.nw import global_align
from scripts.sw import local_align
from scripts.swBlosum import local_align_blosum
from scripts.nwBlosum import global_align_blosum

parser = argparse.ArgumentParser()

parser.add_argument('sequences')
parser.add_argument('method')
parser.add_argument('gap_penalty', nargs='?', default=-1)
parser.add_argument('mismatch_score', nargs='?', default=-1)
parser.add_argument('match_score', nargs='?', default=2)


parser.add_argument('-s','--save', action="store_true")
parser.add_argument('-n','--no_output', action='store_true')
parser.add_argument('-b45','--blosum45', action='store_true')
parser.add_argument('-b50','--blosum50', action='store_true')
parser.add_argument('-b62','--blosum62', action='store_true')
parser.add_argument('-b80','--blosum80', action='store_true')
parser.add_argument('-b90','--blosum90', action='store_true')
args = parser.parse_args()

method = args.method
sequences = Path(args.sequences)
gap_penalty = int(args.gap_penalty)
mismatch_score = int(args.mismatch_score)
match_score = int(args.match_score)
store_results = args.save
no_print_results = args.no_output
b45 = args.blosum45
b50 = args.blosum50
b62 = args.blosum62
b80 = args.blosum80
b90 = args.blosum90


if not sequences.exists():
    print("The target file does not exist")
    raise SystemExit(1)

if b45 == False and b50 == False and b62 == False and b80 == False and b90 == False:
    if method == 'global':
        seq1, seq2, matchlines, seq1_ID, seq2_ID, identity, positives, gaps = global_align(sequences,match_score,mismatch_score,gap_penalty)
        numlines = 70
        row = 0
        len_seq1 = len(seq1)

        if no_print_results == False:

            print(seq1_ID)
            print('|')
            print(seq2_ID)
            print(f'identity: {identity}/{len_seq1}, positives: {positives}/{len_seq1}, gaps: {gaps}/{len_seq1}')

            for i in range(len(seq1)):
                if i % numlines == 0:
                    if (row+1)*numlines < len_seq1:
                        print(row*numlines,(row+1)*numlines)
                    else:
                        print(row*numlines,len_seq1)
                    print(seq1[(row*numlines):((row+1)*numlines)])
                    print(matchlines[(row*numlines):((row+1)*numlines)])
                    print(seq2[(row*numlines):((row+1)*numlines)])
                    row +=1
                i+=1
        
        if store_results == True:
            f = open('output.txt', 'w')
            f.write(f'Alignment of sequences {seq1_ID} and {seq2_ID}\nMethod: Global\nMatch score: {match_score}, mismatch score: {mismatch_score}, gap penalty: {gap_penalty}\n')
            f.write(f'identity: {identity}/{len_seq1}, positives: {positives}/{len_seq1}, gaps: {gaps}/{len_seq1}\n')
            f.close
            row = 0
            for i in range(len(seq1)):
                if i % numlines == 0:
                    f= open('output.txt', 'a')
                    if (row+1)*numlines < len_seq1:
                        f.write(str(row*numlines)+ '    ' + str((row+1)*numlines) +'\n')
                    else:
                        f.write(str(row*numlines)+ '    ' + str(len_seq1) +'\n')
                    f.write(seq1[(row*numlines):((row+1)*numlines)]+'\n')
                    f.write(matchlines[(row*numlines):((row+1)*numlines)]+'\n')
                    f.write(seq2[(row*numlines):((row+1)*numlines)]+'\n')
                    f.close
                    row +=1
                i+=1
            print('Allignment successfull')
            print('Output saved to file \'output.txt\'')

    elif method == 'local':
        seq1, seq2, matchlines, seq1_ID, seq2_ID, start_seq1, start_seq2, identity, positives, gaps = local_align(sequences,match_score,mismatch_score,gap_penalty)
        len_seq1 = len(seq1)
        numlines = 70
        row = 0
        if no_print_results == False:

            print(seq1_ID)
            print('|')
            print(seq2_ID)
            print(f'identity: {identity}/{len_seq1}, positives: {positives}/{len_seq1}, gaps: {gaps}/{len_seq1}')

            for i in range(len(seq1)):
                if i % numlines == 0:
                    if (row+1)*numlines < len_seq1:
                        print(row*numlines + start_seq1,(row+1)*numlines + start_seq1)
                    else:
                        print(row*numlines + start_seq1 ,len_seq1 + start_seq1)
                    print(seq1[(row*numlines):((row+1)*numlines)])
                    print(matchlines[(row*numlines):((row+1)*numlines)])
                    print(seq2[(row*numlines):((row+1)*numlines)])
                    if (row+1)*numlines < len_seq1:
                        print(row*numlines + start_seq2,(row+1)*numlines + start_seq2)
                    else:
                        print(row*numlines + start_seq2 ,len_seq1 + start_seq2)
                    row +=1
                i+=1
        if store_results == True:
            f = open('output.txt', 'w')
            f.write(f'Alignment of sequences {seq1_ID} and {seq2_ID}\nMethod: Local\nMatch score: {match_score}, mismatch score: {mismatch_score}, gap penalty: {gap_penalty}\n')
            f.write(f'identity: {identity}/{len_seq1}, positives: {positives}/{len_seq1}, gaps: {gaps}/{len_seq1}\n')
            f.close
            row = 0

            for i in range(len(seq1)):
                if i % numlines == 0:
                    f= open('output.txt', 'a')
                    if (row+1)*numlines < len_seq1:
                        f.write(str(row*numlines +start_seq1)+ '    ' + str((row+1)*numlines+start_seq1) +'\n')
                    else:
                        f.write(str(row*numlines+start_seq1)+ '    ' + str(len_seq1+start_seq1) +'\n')
                    f.write(seq1[(row*numlines):((row+1)*numlines)]+'\n')
                    f.write(matchlines[(row*numlines):((row+1)*numlines)]+'\n')
                    f.write(seq2[(row*numlines):((row+1)*numlines)]+'\n')
                    if (row+1)*numlines < len_seq1:
                        f.write(str(row*numlines +start_seq2)+ '    ' + str((row+1)*numlines+start_seq2) +'\n')
                    else:
                        f.write(str(row*numlines+start_seq2)+ '    ' + str(len_seq1+start_seq2) +'\n')
                    f.close
                    row +=1
                i+=1
        print('Allignment successfull')
        if store_results == True:
            print('Output saved to file \'output.txt\'')
    else:
        print("Method does not exist.", method)
        raise SystemExit(1)
    
else:
    matrix = 0
    if b45 == True:
        matrix = 45
    elif b50 == True:
        matrix = 50
    elif b62 == True:
        matrix = 62
    elif b80 == True:
        matrix = 80
    elif b90 == True:
        matrix = 90
    

    if method == 'global':
        seq1, seq2, matchlines, seq1_ID, seq2_ID, identity, positives, gaps = global_align_blosum(sequences, matrix)
        numlines = 70
        row = 0
        len_seq1 = len(seq1)

        if no_print_results == False:

            print(seq1_ID)
            print('|')
            print(seq2_ID)
            print(f'identity: {identity}/{len_seq1}, positives: {positives}/{len_seq1}, gaps: {gaps}/{len_seq1}')

            for i in range(len(seq1)):
                if i % numlines == 0:
                    if (row+1)*numlines < len_seq1:
                        print(row*numlines,(row+1)*numlines)
                    else:
                        print(row*numlines,len_seq1)
                    print(seq1[(row*numlines):((row+1)*numlines)])
                    print(matchlines[(row*numlines):((row+1)*numlines)])
                    print(seq2[(row*numlines):((row+1)*numlines)])
                    row +=1
                i+=1
            
        if store_results == True:
            f = open('output.txt', 'w')
            f.write(f'Alignment of sequences {seq1_ID} and {seq2_ID}\nMethod: Global\nSubstitution matrix: BLOSUM{matrix}\n')
            f.write(f'identity: {identity}/{len_seq1}, positives: {positives}/{len_seq1}, gaps: {gaps}/{len_seq1}\n')
            f.close
            row = 0
            for i in range(len(seq1)):
                if i % numlines == 0:
                    f= open('output.txt', 'a')
                    if (row+1)*numlines < len_seq1:
                        f.write(str(row*numlines)+ '    ' + str((row+1)*numlines) +'\n')
                    else:
                        f.write(str(row*numlines)+ '    ' + str(len_seq1) +'\n')
                    f.write(seq1[(row*numlines):((row+1)*numlines)]+'\n')
                    f.write(matchlines[(row*numlines):((row+1)*numlines)]+'\n')
                    f.write(seq2[(row*numlines):((row+1)*numlines)]+'\n')
                    f.close
                    row +=1
                i+=1
            print('Allignment successfull')
            print('Output saved to file \'output.txt\'')

    elif method == 'local':
        seq1, seq2, matchlines, seq1_ID, seq2_ID, start_seq1, start_seq2, identity, positives, gaps = local_align_blosum(sequences,matrix)
        len_seq1 = len(seq1)
        numlines = 70
        row = 0
        if no_print_results == False:

            print(seq1_ID)
            print('|')
            print(seq2_ID)
            print(f'identity: {identity}/{len_seq1}, positives: {positives}/{len_seq1}, gaps: {gaps}/{len_seq1}')

            for i in range(len(seq1)):
                if i % numlines == 0:
                    if (row+1)*numlines < len_seq1:
                        print(row*numlines + start_seq1,(row+1)*numlines + start_seq1)
                    else:
                        print(row*numlines+ start_seq1,len_seq1 + start_seq1)
                    print(seq1[(row*numlines):((row+1)*numlines)])
                    print(matchlines[(row*numlines):((row+1)*numlines)])
                    print(seq2[(row*numlines):((row+1)*numlines)])
                    if (row+1)*numlines < len_seq1:
                        print(row*numlines + start_seq2,(row+1)*numlines + start_seq2)
                    else:
                        print(row*numlines+ start_seq2,len_seq1 + start_seq2)
                    print('--------')
                    row +=1
                i+=1
        if store_results == True:
            f = open('output.txt', 'w')
            f.write(f'Alignment of sequences {seq1_ID} and {seq2_ID}\nMethod: Local\nSubstitution matrix: BLOSUM{matrix}\n')
            f.write(f'identity: {identity}/{len_seq1}, positives: {positives}/{len_seq1}, gaps: {gaps}/{len_seq1}\n')
            f.close
            row = 0

            for i in range(len(seq1)):
                if i % numlines == 0:
                    f= open('output.txt', 'a')
                    if (row+1)*numlines < len_seq1:
                        f.write(str(row*numlines +start_seq1)+ '    ' + str((row+1)*numlines+start_seq1) +'\n')
                    else:
                        f.write(str(row*numlines+start_seq1)+ '    ' + str(len_seq1+start_seq1) +'\n')
                    f.write(seq1[(row*numlines):((row+1)*numlines)]+'\n')
                    f.write(matchlines[(row*numlines):((row+1)*numlines)]+'\n')
                    f.write(seq2[(row*numlines):((row+1)*numlines)]+'\n')
                    if (row+1)*numlines < len_seq1:
                        f.write(str(row*numlines +start_seq2)+ '    ' + str((row+1)*numlines+start_seq2) +'\n')
                    else:
                        f.write(str(row*numlines+start_seq2)+ '    ' + str(len_seq1+start_seq2) +'\n')
                    f.close
                    row +=1
                i+=1
        print('Allignment successfull')
        if store_results == True:
            print('Output saved to file \'output.txt\'')
    else:
        print("Method does not exist.", method)
        raise SystemExit(1)
