# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 10:41:53 2023

@author: Administrator
"""
import numpy as np
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

plt.figure(dpi=600)
plt.scatter(hifdot, alldot, s=3)
plt.xlabel('dE_HF', fontsize=16)
plt.ylabel('dE_T', fontsize=16)

def a_given_point(x, y, cutoff):
    dis = [np.sqrt(8**2*(hifdot[i]-x)**2 + (alldot[i]-y)**2)
           for i in range(len(alldot))]
    sorted_indices = np.argsort(dis)
    idx = sorted_indices[:cutoff]
    return(idx)


def a_given_point_dis(x, y, cutoff):
    dis = [np.sqrt(8**2*(hifdot[i]-x)**2 + (alldot[i]-y)**2)
           for i in range(len(alldot))]
    return(min(dis))


