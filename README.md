# AFF2024

```python
#   █████╗ ███████╗███████╗██████╗  ██████╗ ██████╗ ██╗  ██╗
#  ██╔══██╗██╔════╝██╔════╝╚════██╗██╔═████╗╚════██╗██║  ██║
#  ███████║█████╗  █████╗   █████╔╝██║██╔██║ █████╔╝███████║
#  ██╔══██║██╔══╝  ██╔══╝  ██╔═══╝ ████╔╝██║██╔═══╝ ╚════██║
#  ██║  ██║██║     ██║     ███████╗╚██████╔╝███████╗     ██║
#  ╚═╝  ╚═╝╚═╝     ╚═╝     ╚══════╝ ╚═════╝ ╚══════╝     ╚═╝
#
#  ┏┓┓  ┓   ┏┓  ┓ ┓┏┓  ╻  ┏┓           •
#  ┣┫┃┏┓┣┓┏┓┣ ┏┓┃┏┫┏┛ ━╋━ ┣ ┏┓┓┏┏╋┏┓┏┓╋┓┏┓┏┓
#  ┛┗┗┣┛┛┗┗┻┻ ┗┛┗┗┻┗━  ╹  ┻ ┛ ┗┛┛┗┛ ┗┻┗┗┗┛┛┗
#     ┛
#       Predicting protein conformational motions
#  using energetic frustration analysis and AlphaFold2
```



## Introduction

This Github repository is for our work: *Predicting protein conformational motions using energetic frustration analysis and AlphaFold2*

We have prepared folders containing the scripts and data in this work. The structures and descriptions are listed as following.

## Walkthrough

We also provide a notebook for a quick walkthrough: `Walkthrough.ipynb`[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Gxy-with-luv/AFF2024/blob/main/Walkthrough.ipynb), with this one can get a frustration-energy landscape mapped from the MSA space using physically-based frustration information for AdK (with closed conformation as the native state) and then generate the alternative open state using High $\Delta E_{HF}$ region, masking and mixing methods.

## Folder structures and descriptions 

### Core prediction part:

Using Colab is a convenient way to install and load AF2-related modules without the burden of configurations for environment. We provide a notebook: `Main_Prediction.ipynb` to load and use AlphaFold2, to predict structures from the MSA hints given by ourselves.

### AdK Folder:

``` bash
AdK # Every sub-folder also contain scripts for analysis and output results
├── Frustration-Energy-space # Mapping the MSA space onto the Frustration-Energy surface
│   ├── circ50 # Sliding from low dE_{HF} region towards high dE_{HF} region
│   ├── circ50higher # Sliding upwards to explore the energy surface learnt by AF2
│   ├── sampling400 # Random sampling for structure predictions
│   └── ...
│
├── AdK-erasing # Erasing/Masking highly-frustrated sites pushes the structure towards the open state
│   ├── erase_highFrus # top 30 with ratio 0.5, masking the most highly-frustrated sites
│   ├── erase_highFrus_various_ratio # Masking the most highly-frustrated sites with various masking ratio
│   ├── erase_lowFrus # Minimally-frustrated sites, as a control for comparison.
│   ├── erase_lowFrus_various_ratio # Minimally-frustrated sites, as a control for comparison.
│   ├── erase_nobreakage # Sites with no breakage during conformational change, as a control for comparison.
│   ├── erase_breakage # Sites with breakage during conformational change, as a control for comparison.
│   └── ...
│
├── AdK-mixing # Mixing for generating the conformational-changing pathway
│   ├── Mixing100 # Mixing with size of 100 sequences for the Open-state-set
│   ├── Mixing50 # Mixing with size of 50 sequences for the Open-state-set
│   └── ...
│
├── AdK-PCA_analysis # Containing PCA analysis for coordinates, distances(strain) and contacts(cracking) respectively 
│
└── SI-AA_simulation # Gromacs setup files and also two output trajectories

```

###  Rosetta scripts Folder:

```bash
*Rosetta Folder:*

Rosetta_scripts
├── Rosetta_AdK
├── Rosetta_KaiB
└── Rosetta_RBP

*Under each script folder:*

Rosetta_AdK
├── native.pdb # Native structure
├── native.seq # Query sequence
├── native.xml # Rosetta scripts to calculate the native energy
├── test.xml # Rosetta scripts to calculate decoy energies
├── *.sh # Calculation job submission
├── peps # Containing sequences from the MSA space
├── *.py # Scripts for collecting the energy profiles from Rosetta calculation
└── *.npy # Output energy profiles
```



