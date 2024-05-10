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

def traj2con(traj, length):
    distances, pairs = md.compute_contacts(traj, 
                                           contacts='all', 
                                           scheme='closest-heavy', 
                                           ignore_nonprotein=False)
    threshold = 1.0
    
    # Find contacts within the threshold
    contacts_within_threshold = distances < threshold
    
    # Filter pairs based on the threshold
    contact_pairs = pairs[contacts_within_threshold.flatten()]
    
    contacts = np.zeros([length, length])
    for pair in contact_pairs:
        contacts[pair[0], pair[1]] = 1
        contacts[pair[1], pair[0]] = 1
    return(contacts)

traj = md.load('ensemble.pdb')
traj = traj.atom_slice(traj.top.select('name CA'))
contacts_colle = []
for frame in traj:
    contacts = traj2con(frame, length=214)
    contacts_colle.append(contacts)

combined_distances = np.array(contacts_colle).reshape([-1, 214*214])

from sklearn.decomposition import PCA

pca = PCA(n_components=20)  
pca.fit(combined_distances)

transformed_distances = pca.transform(combined_distances)

explained_variance_ratio = pca.explained_variance_ratio_

principal_components = pca.components_

a = principal_components[0] #change here
a = a.reshape([214,214])

