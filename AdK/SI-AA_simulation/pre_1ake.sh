#!/bin/env bash
module purge
module load compiler/cmake-compiler/3.20.1
module load compiler/gcc-compiler/7.3.1
module load compiler/intel-compiler/2021.3.0
module load compiler/cuda/11.4-with-cuDNN8.2.4

source /public/home/win1403/software/gromacs_gpu_plumed/bin/gmx-completion.bash
source /public/home/win1403/software/gromacs_gpu_plumed/bin/gmx-completion-gmx.bash
source /public/home/win1403/software/gromacs_gpu_plumed/bin/GMXRC
export MODULEPATH=/public/home/win1403/software/plumed_gpu/lib/plumed:$MODULEPATH   # cpu + plumed
module load modulefile

pdb_name="1ake_A"
topology_file="topology.top"
mdp_ions="/public/home/win1403/data/AF_cluster/ADK/pre_file/ions.mdp"
force_field='amber99sb-ildn'
ions_concentration=0.15

# choose force-field
gmx pdb2gmx -f ${pdb_name}.pdb -o ${pdb_name}.gro -ff ${force_field} -water spce -p ${topology_file}

# Solvate the protein
gmx editconf -f ${pdb_name}.gro -o boxed.gro -c -d 1.0 -bt cubic

gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p ${topology_file}

# Neutralize the system
gmx grompp -f ${mdp_ions} -c solvated.gro -p ${topology_file} -o ions.tpr
echo 13 | gmx genion -s ions.tpr -o solvated_ions.gro -p ${topology_file} -pname NA -nname CL -neutral -conc ${ions_concentration}

