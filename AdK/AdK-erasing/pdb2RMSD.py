#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 10:23:22 2024

@author: tenran
"""

import numpy as np
import mdtraj as md

dirname = 'erase_0contact_easy_results/'

import glob

#pdbs = glob.glob(dirname+'harder*pdb')
pdbs = glob.glob(dirname+'*pdb')
closedca = md.load_pdb('../closeCA.pdb')
openca = md.load_pdb('../openCA.pdb')

box = []
for pdb in pdbs:
    pdb = md.load_pdb(pdb)
    pdb = pdb.atom_slice(pdb.top.select('name CA'))
    rmsdc = md.rmsd(pdb, closedca)[0]
    rmsdo = md.rmsd(pdb, openca)[0]
    box.append(rmsdc)

import matplotlib.pyplot as plt
plt.violinplot(box)    