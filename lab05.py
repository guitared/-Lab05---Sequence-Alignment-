import argparse,itertools
import sys,os
def print_matrix(matrix,seq1,seq2):
    seq2 = '  '+seq2
    matrix.insert(0,list(' '+seq1))
    matrix = [[seq2[i]] + x for i,x in enumerate(matrix)]
    print
    for row in matrix:
        if row[0]=='-': continue
        for i,val in enumerate(row):
            if matrix[0][i]=='-': break
            print '{:>3}'.format(val),
        print
        print
def print_backtrack(matrix,seq1,seq2):
    matrix[0][0]=' X '
    seq2 = '  '+seq2
    matrix.insert(0,list(' '+seq1))
    matrix = [[seq2[i]] + x for i,x in enumerate(matrix)]
    print
    for row in matrix:
        if row[0]=='-': continue
        for i,val in enumerate(row):
            if matrix[0][i]=='-': break
            if val==1: val='  |'
            if val==2: val='_  '
            if val==3: val='_ |'
            if val==4: val=' \\ '
            if val==5: val=' \\|'
            if val==6: val='_\\ '
            if val==7: val='_\\|'
            print '{:>3}'.format(val),
        print
        print
def global_backtrack(matrix,backtrackmatrix,seq1,seq2,b1,b2):
    alignments = []
    i=len(seq2)-len(b2)
    j=len(seq1)-len(b1)
    while len(b2)>len(b1):
        b1+='-'
    while len(b1)>len(b2):
        b2+='-'
    c=0
    #print '--------------------------------'
    while j>0 and i>0:
        #print c,i,j
        c+=1
        if backtrackmatrix[i][j] == 0b001:
            #print "del"
            b1+='-'
            b2+=seq2[i-1]
            i-=1
        elif backtrackmatrix[i][j] == 0b010:
            #print "ins"
            b1+=seq1[j-1]
            b2+='-'
            j-=1
        elif backtrackmatrix[i][j] == 0b100:
            #print "mat"
            b1+=seq1[j-1]
            b2+=seq2[i-1]
            i-=1
            j-=1
        elif backtrackmatrix[i][j] == 0b011:
            b3,b4 = '%s' % b1, '%s' % b2
            b1+='-'
            b2+=seq2[i-1]
            b3+=seq1[j-1]
            alignments.append(global_backtrack(matrix,backtrackmatrix,seq1,seq2,b3,b4))
            i-=1
        elif backtrackmatrix[i][j] == 0b101:
            b3,b4 = '%s' % b1, '%s' % b2
            b1+=seq1[j-1]
            b2+=seq2[i-1]
            b4+=seq2[i-1]
            alignments.append(global_backtrack(matrix,backtrackmatrix,seq1,seq2,b3,b4))
            i-=1
            j-=1
        elif backtrackmatrix[i][j] == 0b110:
            b3,b4 = '%s' % b1, '%s' % b2
            b1+=seq1[j-1]
            b2+=seq2[i-1]
            b3+=seq1[j-1]
            alignments.append(global_backtrack(matrix,backtrackmatrix,seq1,seq2,b3,b4))
            i-=1
            j-=1
        elif backtrackmatrix[i][j] == 0b111:
            b3,b4 = '%s' % b1, '%s' % b2
            b5,b6 = '%s' % b1, '%s' % b2
            b1+=seq1[j-1]
            b2+=seq2[i-1]
            b3+=seq1[j-1]
            b6+=seq2[i-1]
            alignments.append(global_backtrack(matrix,backtrackmatrix,seq1,seq2,b3,b4))
            alignments.append(global_backtrack(matrix,backtrackmatrix,seq1,seq2,b5,b6))
            i-=1
            j-=1
        #print "b1="+b1[::-1]
        #print "b2="+b2[::-1]
    i=len(seq2)-len(b2)
    while i>0:
        b2+=seq2[i-1]
        b1+='-'
        i-=1
    i=len(seq1)-len(b1)
    while i>0:
        b1+=seq1[i-1]
        b2+='-'
        i-=1
    alignments.append((b1[::-1],b2[::-1]))
    #print '--------------------------------'
    return list(itertools.chain.from_iterable(alignments))
def global_alignment(seq1,seq2,match,mismatch,gap,show_score,show_backtrack):
    print 'Global alignment'
    print 'sequence#1: '+seq1
    print 'sequence#2: '+seq2
    matrix = [[0 for y in range(len(seq1)+1)] for x in range(len(seq2)+1)]
    backtrackmatrix = [[0 for y in range(len(seq1)+1)] for x in range(len(seq2)+1)]
    for i in range(len(seq1)):
        matrix[0][i+1]=matrix[0][i]-1
        backtrackmatrix[0][i+1]=0b010
    for i in range(len(seq2)):
        matrix[i+1][0]=matrix[i][0]-1
        backtrackmatrix[i+1][0]=0b001
    for j in range(len(seq1)):
        for i in range(len(seq2)):
            if seq1[j]==seq2[i]:
                matc = matrix[i][j] + match
            elif seq1[j]=='-' or seq2[i]=='-':
                matc = matrix[i][j] + gap
            else:
                matc = matrix[i][j] + mismatch
            inst = matrix[i+1][j] + gap
            dele = matrix[i][j+1] + gap
            maxx = max(matc,inst,dele)
            if matc == maxx:
                backtrackmatrix[i+1][j+1]+=0b100
            if inst == maxx:
                backtrackmatrix[i+1][j+1]+=0b010
            if dele == maxx:
                backtrackmatrix[i+1][j+1]+=0b001
            matrix[i+1][j+1]= maxx
    alignments = global_backtrack(matrix,backtrackmatrix,seq1,seq2,'','')
    if show_score: print_matrix(matrix,seq1,seq2)
    if show_backtrack: print_backtrack(backtrackmatrix,seq1,seq2)
    print "\nBest alignment(s)\n"
    for i,align in enumerate(alignments):
        print align
        if i%2==1:
            print

def local_alignment(seq1,seq2,match,mismatch,gap):
    print 'local alignment'
    print 'sequence#1: '+seq1
    print 'sequence#2: '+seq2
    matrix = [[0 for x in range(len(seq1)+1)] for x in range(len(seq2)+1)]
    print_matrix(matrix,seq1,seq2)


parser = argparse.ArgumentParser(description='Lab05 - Sequence Alignment using dynamic programming Assignment')
#input argument
input_group = parser.add_mutually_exclusive_group(required=True)
input_group.add_argument('-i','--input',metavar='name',help='input FASTA file name',type=argparse.FileType('r'),default=sys.stdin)
input_group.add_argument('-S','--sequnces',metavar='seq',nargs=2,help='input two sequences')
# options argument
parser.add_argument('-s','--scoringmatrix',help='show scoring matrix',action="store_true")
parser.add_argument('-b','--backtrackmatrix',help='show backtracking pointer matrix',action="store_true")
# alignment parameter argument
parser.add_argument('-m','--match',metavar='num',help='match value (default:1)',type=int,default=1)
parser.add_argument('-M','--mismatch',metavar='num',help='mismatch value (default:-1)',type=int,default=-1)
parser.add_argument('-g','--gap',metavar='num',help='gap value (default:-1)',type=int,default=-1)
# output argument
parser.add_argument('-o','--output',metavar='name',default=False,help='generate output to file with name')
args = parser.parse_args()

# processes input data
# ------------------------------------------------------------------------------
input_seq = ['','']
if args.sequnces:
    input_seq=args.sequnces
elif args.input:
    data = args.input.read()
    i=-1
    for line in data.split('\n'):
        if line.startswith(';'):
            continue
        if line.startswith('>'):
            i+=1
            if i==2: # get only first two sequcne in file
                break
        else:
            input_seq[i]+=line
global_alignment(input_seq[0],input_seq[1],args.match,args.mismatch,args.gap,args.scoringmatrix,args.backtrackmatrix)

if args.output:
    f = open(os.path.splitext(args.output)[0]+'.fas', 'w')
