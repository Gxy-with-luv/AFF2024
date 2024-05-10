#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:59:56 2024

@author: tenran
"""

import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt


def a3m_seq_keep(msafile):
    seqs = []
    with open(msafile, 'r') as fi:
        for line in fi:
            if '>' not in line:
                seqs.append(line)
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


allseqs = a3m_seq_keep('../MSA.a3m')

scores = np.load('../ENE_sing.npy')
scores = scores.reshape([-1, 214])
natscore = np.load('../ENE_nat.npy')
natscore = natscore.reshape([-1, 214])

scores = scores - natscore

filename = '../frustration_configurational_5adens'
frustbox = []
with open(filename, 'r') as fi:
    for line in fi:
        if (len(line) > 5) and ('Residue' not in line):
            line = line.split()
            line = np.array(line[-3:]).astype(np.float32)
            frustbox.append(line)

frustbox = np.array(frustbox)

hif = frustbox[:, 0]
lof = frustbox[:, 1]
hifdot = np.dot(scores, hif)
alldot = np.sum(scores, axis=1)

plt.scatter(hifdot, alldot)


boxxy = []
rmsd = []
closedpdb = md.load_pdb('../closedCA.pdb')
openpdb = md.load_pdb('../openCA.pdb')

for roll in range(0,400):
    a3m = '../sampling400/samp'+str(roll)+'.a3m'
    seqs = a3m_seq_keep(a3m)
    idxs = [allseqs.index(i) for i in seqs]
    x = np.mean(hifdot[idxs])
    y = np.mean(alldot[idxs])
    boxxy.append([x,y])
    pdb = '../Predicted_results/sampling400_results/samp'+str(roll)+'.pdb'
    pdb = md.load_pdb(pdb)
    pdb = pdb.atom_slice(pdb.top.select('name CA'))
    rmsdc = md.rmsd(pdb, closedpdb)[0]
    rmsdo = md.rmsd(pdb, openpdb)[0]
    rmsd.append([rmsdc, rmsdo])

boxxy = np.array(boxxy)
rmsd = np.array(rmsd)

plt.figure(dpi=600)
plt.scatter(boxxy[:, 0], boxxy[:, 1], c=rmsd[:, 1])
plt.colorbar()
plt.xlabel('HiF_Energy')
plt.ylabel('Total_Energy')
