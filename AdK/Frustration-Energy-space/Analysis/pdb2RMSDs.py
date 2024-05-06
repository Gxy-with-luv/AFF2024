# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 11:10:44 2023

@author: Administrator
"""
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

allbox1 = []
allbox2 = []
allpdb = []
hist = []

# Please replace the pdbdir and numbers with corresponding name and value.
pdbdir = '../Predicted_results/circ50_results/'
bigbox = []
openpdb = md.load_pdb('../openCA.pdb')
closedpdb = md.load_pdb('../closedCA.pdb')
for roll in [0]:
    box1 = []
    box2 = []
    for num in np.arange(60):
        num = str(num)
        pdbfile = pdbdir + 'circ'+str(num)+ '.pdb'
        pdb = md.load_pdb(pdbfile)
        allpdb.append(pdb)
        atoms = pdb.top.select('name CA')
        pdb = pdb.atom_slice(atoms)
        pdb.superpose(openpdb, 0)
        d1 = md.rmsd(pdb, openpdb)
        pdb.superpose(closedpdb, 0)
        d2 = md.rmsd(pdb, closedpdb)
        box1.append(d1[0])
        box2.append(d2[0])
    allbox1.extend(box1)
    allbox2.extend(box2)

