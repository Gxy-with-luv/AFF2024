# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 21:27:08 2023

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


natseq = "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"
allseqs = a3m_seq_keep('MSA.a3m')

scores = np.load('ENE_sing.npy')
scores = scores.reshape([-1, 214])
natscore = np.load('ENE_nat.npy')
natscore = natscore.reshape([-1, 214])

scores = scores - natscore

filename = 'frustration_configurational_5adens'
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

hifdot = np.dot(scores, hif)
lofdot = np.dot(scores, lof)
alldot = np.sum(scores, axis=1)


def a_given_point(x, y):
    plt.figure(dpi=300)
    dis = [np.sqrt(7.7**2*(hifdot[i]-x)**2 + (alldot[i]-y)**2)
           for i in range(len(alldot))]
    sorted_indices = np.argsort(dis)
    idx = sorted_indices[:50]
    plt.scatter(hifdot, alldot, s=2)
    plt.xlabel('HiFrustration Score')
    plt.ylabel('Energy Score')
    plt.scatter(hifdot[idx], alldot[idx], c='r', s=2)
    plt.show()
    return(idx)


def a_given_point_min(x, y):
    # To be used to see whether a sampling is successful in [Random sampling 400 times]
    dis = [np.sqrt(7.7**2*(hifdot[i]-x)**2 + (alldot[i]-y)**2)
           for i in range(len(alldot))]
    sorted_indices = np.argsort(dis)
    return(min(dis))


# Grep the frustration-energy space (dEhf_dEt)
# This is Fig2C
# #############################################
plt.figure(dpi=600)
plt.scatter(hifdot, alldot, s=3)
plt.xlabel('dEhf', fontsize=16)
plt.ylabel('dEt', fontsize=16)
plt.show()
# #############################################


# Sliding Results (dEhf_dEt)
# This is Fig3C
# #############################################
for i in range(60):
    num = i/2-5
    idx = a_given_point(num, 0)
    with open('circ50/circ'+str(i)+'.a3m', 'w') as fo:
        fo.write('>101\n')
        fo.write("MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG\n")
        count = 0
        for j in idx:
            fo.write('>'+str(count)+'\n')
            fo.write(allseqs[j])
# #############################################


# Random sampling Results upon Higher dEt
# This is Fig5A
# #############################################
for ene_all in [20, 30, 40, 50, 60, 70, 80, 90, 100, 110]:
    idxbox = [i for i in range(len(alldot)) if np.abs(
        alldot[i]-ene_all) < 10 and hifdot[i] > 0]
    plt.scatter(hifdot, alldot)
    for i in range(30):
        idx = np.random.choice(idxbox, size=50, replace=False)
        with open('circ50higher/circ'+str(ene_all)+'_'+str(i)+'.a3m', 'w') as fo:
            fo.write('>101\n')
            fo.write("MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG\n")
            count = 0
            for j in idx:
                fo.write('>'+str(count)+'\n')
                fo.write(allseqs[j])
                count += 1
# #############################################


# Random sampling 400 times
# This is Fig2E
# #############################################
boxxy = []
for i in range(400):
    ifsuccess = False
    while ifsuccess == False:
        x = np.random.random()*27-7
        y = np.random.random()*150-20
        if a_given_point_min(x, y) < 2:
            ifsuccess = True
            boxxy.append([x, y])
        idxs = a_given_point(x, y)
        with open('sampling400/samp'+str(i)+'.a3m', 'w') as fo:
            fo.write('>101\n')
            fo.write("MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG\n")
            count = 0
            for j in idxs:
                fo.write('>'+str(count)+'\n')
                fo.write(allseqs[j])
                count += 1
# #############################################
