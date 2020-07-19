import sys
import os
import pandas as pd
from Bio import SeqIO

def select(i, B):
    Bl = list(str(B))
    cnt = 0
    for p,x in enumerate(Bl):
        if x == '1':
            cnt += 1
            if cnt == i:
                return p

def rank(i, B):
    Bl = list(str(B))
    return Bl[:i].count('1')

def get_S(j, B, I):
    rnk_ind = rank(j, I)
    print("rank I=",rnk_ind)
    j_zero_based = j - 1
    if list(I)[j_zero_based] == '1':
        print('inside')
        return(M[rnk_ind-1])
    rnk = rank(j, B)
    print("rank B=",rnk)
    sel = select(rnk, B)
    print("select=", sel)
    init = Q[rnk-1] - 1
    print("initial pos of phrase=", init)
    cs = j - sel - 1
    print(init+cs)
    if (j % 1000) == 0:
        print(j)
    return R[init + cs]


if __name__ == "__main__":

    #R = sys.argv[1]
    #S   = sys.argv[2]
    #M = 'T C T C A C C C T T C C T G C T G A A T C G A A C A T A T G A A A T'.split(" ")
    for reca in SeqIO.parse(sys.argv[1], 'fasta'):
        R = str(reca.seq)

    for recb in SeqIO.parse(sys.argv[2], 'fasta'):
        S = str(recb.seq)

    # with open(sys.argv[1], 'r') as infile:
    #     for l in infile.readlines():
    #         R = l
    #
    # with open(sys.argv[2], 'r') as infile:
    #     for l in infile.readlines():
    #         S = l

    with open('outM.txt', 'r') as outMf:
        for l in outMf.readlines():
            M = l.split(" ")

    with open('outI.txt', 'r') as outIf:
        for l in outIf.readlines():
            B = l[:-1]
    with open('outB.txt', 'r') as outBf:
        for l in outBf.readlines():
            I = l[:-1]


    df = pd.read_csv('out_sdsl.txt', sep="\t", header=None)
    df = df.drop(df.index[0])
    df = df.loc[~(df.iloc[:,1] == '$'), :]
    Q = [int(x) for x in pd.to_numeric(df.iloc[:, 1], downcast='integer').values + 1]


    reconstructed_S = [get_S(x+1, B, I) for x in range(len(S))]
    recon_S = ''.join(reconstructed_S)

    if(recon_S == S):
        print("SUCCESS !")
    else:
        print("ERROR")
