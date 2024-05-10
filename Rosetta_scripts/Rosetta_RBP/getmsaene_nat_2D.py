# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:01:04 2023

@author: Administrator
"""
import numpy as np

pairs = []
for i in range(0, 271):
    for j in range(i, 271):
        pairs.append([i, j])
box = []
enes = []
singres = []
enelog = []
logfile = 'native.log'
eneres = np.zeros([271])
ene2D = np.zeros([271, 271])
with open(logfile, 'r') as fi:
    for line in fi:
        line = line.split()
        if 'Res1' not in line and 'nonzero' not in line and float(line[4]) <= 5.000:
            aa1 = int(float(line[1].split('_')[1][1:])) - 1
            aa2 = int(float(line[2].split('_')[1][1:])) - 1
            ene = float(line[-1]) - 1.0*(float(line[4])) - \
                float(line[10]) - float(line[15])
            eneres[aa1] += 0.5 * ene
            eneres[aa2] += 0.5 * ene
            ene2D[aa1, aa2] = ene
            ene2D[aa2, aa1] = ene
for pair in pairs:
    enelog.append(eneres[pair[0]] + eneres[pair[1]])
enes.append(enelog)
singres.append(eneres)
enes = np.array(enes)
meanene = np.mean(enes, axis=0)
stdene = np.std(enes, axis=0)
singres = np.array(singres)
meansing = np.mean(singres, axis=0)
stdsing = np.std(singres, axis=0)
box.append(singres)
#np.save('natbox.npy', box)
np.save('natbox2D.npy', ene2D)
