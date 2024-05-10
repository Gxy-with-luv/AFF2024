# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 14:55:44 2023

@author: Administrator
"""
import numpy as np

from tqdm import tqdm


def a3m_seq_keep(msafile):
    natseq = "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"  # Example sequence
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

# sorted_indices = np.argsort(hif)[::-1]

contacts = np.load('num_contact.npy')
msafile = 'original50.a3m'
contact0 = [i for i in range(214) if contacts[i]>=3]
contact2 = [i for i in range(214) if contacts[i]==2]
for ratio in [0.5]:
    for roll in np.arange(80):
        lst_temp = np.random.choice(contact2, 6)
        lst = contact0.copy()
        lst.extend(lst_temp)
        ratio = round(ratio, 1)
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
        with open('erase_2contact/erase'+str(ratio)+'_'+str(roll)+'.a3m', 'w') as fo:
            fo.write('>101\n')
            fo.write('MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG\n')
            for seq in range(len(newseqs)):
                fo.write('>'+str(seq)+'\n')
                fo.write(newseqs[seq])
