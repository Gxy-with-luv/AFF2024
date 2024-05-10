#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 10:17:16 2024

@author: tenran
"""

import numpy as np
import mdtraj as md
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

traj = md.load_pdb('ensemble.pdb')
ref_frame = md.load_pdb('closed.pdb')
ref_frame = ref_frame.atom_slice(ref_frame.top.select('name CA'))
traj = traj.atom_slice(traj.top.select('name CA'))
traj = traj.superpose(ref_frame)
coord = traj.xyz.reshape(traj.n_frames, -1)

pca = PCA(n_components=20)
pca.fit(coord)

pca_coords = pca.transform(coord)

plt.scatter(pca_coords[:, 0], pca_coords[:, 1])
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Principal Component Analysis')
plt.show()
transformed_distances = pca.transform(coord)
eigenvectors = pca.components_
tosee = eigenvectors[0].reshape([214,3])

xyz = ref_frame.xyz[0]
with open('pca_coord.txt', 'w') as fi:
    for i in range(214):
        x = xyz[i]
        if np.linalg.norm(tosee[i]) > 0.05:
            xe = 7*tosee[i] + x
            x = ', '.join([str(round(10*j,4)) for j in x])
            xe = ', '.join([str(round(10*j,4)) for j in xe])
            line = 'cgo_arrow ['+x+'], ['+xe+'], radius=0.3, gap=0'
            print(line)
            fi.write(line+'\n')
