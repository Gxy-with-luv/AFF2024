# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 17:34:04 2023

@author: Administrator

Grep the conformational change in AdK
"""

import numpy as np
import mdtraj as md

traj = md.load_pdb('pdbfilename.pdb')

CAs = traj.top.select('name CA')
corelid = [i for i in traj.top.select('resid >= 114') if i in traj.top.select('resid <= 124')]
corelid = [i for i in corelid if i in CAs]
core = [i for i in traj.top.select('resid >= 89') if i in traj.top.select('resid <= 99')]
core = [i for i in core if i in CAs]
nmp = [i for i in traj.top.select('resid >= 34') if i in traj.top.select('resid <= 54')]
nmp = [i for i in nmp if i in CAs]
core2 = [i for i in traj.top.select('resid >= 178') if i in traj.top.select('resid <= 184')]
core2 = [i for i in core2 if i in CAs]
lid = [i for i in traj.top.select('resid >= 124') if i in traj.top.select('resid <= 152')]
lid = [i for i in lid if i in CAs]

corelidxyz = np.mean(np.array([traj.xyz[:,i,:] for i in corelid]), axis=0)
corexyz = np.mean(np.array([traj.xyz[:,i,:] for i in core]), axis=0)
core2xyz = np.mean(np.array([traj.xyz[:,i,:] for i in core2]), axis=0)
nmpxyz = np.mean(np.array([traj.xyz[:,i,:] for i in nmp]), axis=0)
lidxyz = np.mean(np.array([traj.xyz[:,i,:] for i in lid]), axis=0)


ang1s = []
ang2s = []
for frame in range(len(traj)):
    corelid_core = np.linalg.norm(corelidxyz[frame] - corexyz[frame])
    core_nmp = np.linalg.norm(nmpxyz[frame] - corexyz[frame])
    corelid_nmp = np.linalg.norm(corelidxyz[frame] - nmpxyz[frame])
    
    core2_corelid = np.linalg.norm(core2xyz[frame] - corelidxyz[frame])
    corelid_lid = np.linalg.norm(corelidxyz[frame] - lidxyz[frame])
    core2_lid = np.linalg.norm(core2xyz[frame] - lidxyz[frame])
    
    ang1 = np.arccos((corelid_core**2 + core_nmp**2 - corelid_nmp**2)/2/corelid_core/core_nmp)/np.pi*180
    ang2 = np.arccos((core2_corelid**2 + corelid_lid **2 - core2_lid**2)/2/core2_corelid/corelid_lid)/np.pi*180
    
    ang1s.append(ang1)
    ang2s.append(ang2)