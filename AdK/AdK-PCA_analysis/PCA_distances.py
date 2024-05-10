# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 10:24:22 2024

@author: Administrator
"""
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
# Make an ensemble here.

import glob

traj = md.load('ensemble.pdb')
traj = traj.atom_slice(traj.top.select('name CA'))
pairs = [[i,j] for i in range(214)
             for j in range(214)]

distances = md.compute_distances(traj, pairs)

distances = distances.reshape([-1,214,214])

combined_distances = distances.reshape(distances.shape[0], -1)

from sklearn.decomposition import PCA

pca = PCA(n_components=20)  
pca.fit(combined_distances)

transformed_distances = pca.transform(combined_distances)

explained_variance_ratio = pca.explained_variance_ratio_

principal_components = pca.components_

a = principal_components[0] #change here
a = a.reshape([214,214])

native = md.load_pdb('closed.pdb')
native = native.atom_slice(native.top.select('name CA'))
xyz = (native.xyz)[0]
end_box = []

rmsd1s = []
rmsd2s = []
closedCA = md.load('closedCA.pdb')
openCA = md.load('openCA.pdb')
for frame in traj:
    rmsd1 = md.rmsd(frame, closedCA)[0]
    rmsd2 = md.rmsd(frame, openCA)[0]
    rmsd1s.append(rmsd1)
    rmsd2s.append(rmsd2)
    

for i in range(214):
    start = xyz[i]
    end = np.zeros([3])
    for j in range(214):
        if i==j:
            continue
        end_t = -(xyz[j] - xyz[i]) * a[i,j] / np.linalg.norm(xyz[j] - xyz[i])
        end += end_t
    end_box.append(end)
'''
with open('pcatxt2_rev.txt', 'w') as fi:
    for i in range(214):
        x = xyz[i]
        if np.linalg.norm(end_box[i]) > 0.7:
            xe = 4*end_box[i] + x
            x = ', '.join([str(round(10*j,4)) for j in x])
            xe = ', '.join([str(round(10*j,4)) for j in xe])
            line = 'cgo_arrow ['+x+'], ['+xe+'], radius=0.3, gap=0'
            print(line)
            fi.write(line+'\n')'''
    