# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 19:37:23 2023

@author: Administrator
"""

import numpy as np


def a3m_seq_keep(msafile):
    boxx = []
    natseq = "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"  # Example sequence
    seqs = []
    with open(msafile, 'r') as fi:
        for line in fi:
            if '>' not in line:
                seqs.append(line)
    native = seqs[0]
    alter = seqs[1:]
    alter_fin = []
    for alt in alter:
        seq = []
        for i in range(len(alt)):
            if alt[i] == '-':
                seq.append('-')
            else:
                seq.append(alt[i])
        seq = ''.join(seq)
        alter_fin.append(seq)
    return alter_fin


def mix_msa(msafile1, msafile2, ratio):
    seq1 = a3m_seq_keep(msafile1)
    seq2 = a3m_seq_keep(msafile2)
    np.random.shuffle(seq1)
    np.random.shuffle(seq2)
    N1 = len(seq1)
    N2 = len(seq2)
    seq1 = seq1[:int(N1*ratio)]
    seq2 = seq2[:N2 - int(N2*ratio)]
    seq1.extend(seq2)
    np.random.shuffle(seq1)
    return(seq1)


def seq_to_a3m(seq, outfile):
    with open(outfile, 'w') as fo:
        natseq = "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"  # Example sequence
        fo.write('>101\n')
        fo.write(natseq+'\n')
        for i in range(len(seq)):
            fo.write('>'+str(i)+'\n')
            fo.write(seq[i])


msafile1 = 'Closed500.a3m'
msafile2 = 'Open50.a3m'

for ratio in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]:
    for roll in range(10):
        seq = mix_msa(msafile1, msafile2, ratio)
        outfile = 'Mixing100/mix'+str(round(ratio, 3))+'_'+str(roll)+'.a3m'
        seq_to_a3m(seq, outfile)
