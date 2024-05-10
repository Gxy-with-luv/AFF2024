# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 14:55:44 2023

@author: Administrator
"""
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import umap
from sklearn.cluster import KMeans
from Bio.SubsMat import MatrixInfo
from tqdm import tqdm


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

import numpy as np

filename = '../Frustration-Energy-space/frustration_configurational_5adens'
frustbox = []
with open(filename, 'r') as fi:
    for line in fi:
        if (len(line) > 5) and ('Residue' not in line):
            line = line.split()
            line = np.array(line[-3:]).astype(np.float32)
            frustbox.append(line)

frustbox = np.array(frustbox)

hif = frustbox[:, 0]
lof = frustbox[:, -1]

sorted_indices = np.argsort(hif)[::-1]


msafile = 'original50.a3m'

for ratio in [0.5]: # Change this into 0.1~0.9 to reproduce erase_highFrus_various_ratio
    for roll in np.arange(80):
        lst = sorted_indices[:30]
        ratio = round(ratio,1)
        newseqs = []
        seqs = a3m_seq_keep(msafile)
        for seq in seqs:
            dice = np.random.random()
            if dice < ratio:
                seq = list(seq)
                for bar in lst:
                    seq[bar] = '-' 
                seq = ''.join(seq)
            newseqs.append(seq)
        with open('erase_highFrus/erase'+str(ratio)+'_'+str(roll)+'.a3m', 'w') as fo:
            fo.write('>101\n')
            fo.write('MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG\n')
            for seq in range(len(newseqs)):
                fo.write('>'+str(seq)+'\n')
                fo.write(newseqs[seq])
            
