#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 19:29:40 2024

@author: tenran
"""

import numpy as np
import mdtraj as md


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


allseqs = a3m_seq_keep('MSA.a3m')


bigbox = np.load('ENE_sing.npy') 
boxnat = np.load('ENE_nat.npy')
boxnat = boxnat.reshape([-1])
bigbox = bigbox - boxnat

filename = 'KaiBFrus.txt'
frustbox = []
with open(filename, 'r') as fi:
    for line in fi:
        if (len(line) > 5) and ('Residue' not in line):
            line = line.split()
            line = np.array(line[-3:]).astype(np.float32)
            frustbox.append(line)

frustbox = np.array(frustbox)

hif = frustbox[:, 0]
lof = frustbox[:, 2]

hifdot = np.dot(bigbox, hif)
alldot = np.sum(bigbox, axis=1)


box = []
es = md.load_pdb('esCA.pdb')
gs = md.load_pdb('gsCA.pdb')
for i in range(500):
    pdbfile = 'random_sampling_results/circ'+ str(i)+'.pdb'
    pdb = md.load_pdb(pdbfile)
    pdb = pdb.atom_slice(pdb.top.select('name CA'))
    rmsd2 = md.rmsd(pdb, gs)[0]
    rmsd1 = md.rmsd(pdb, es)[0]
    box.append([rmsd2, rmsd1])

boxxy = np.loadtxt('random_sampling_xy.txt')
import matplotlib.pyplot as plt
plt.figure(dpi=600)
plt.scatter(hifdot, alldot, marker='o', s=4, alpha=0.2, c='gray')
plt.scatter(boxxy[:,0], boxxy[:,1], c=np.array(box)[:,0], marker='o', s=10, alpha=0.3)
plt.colorbar()

