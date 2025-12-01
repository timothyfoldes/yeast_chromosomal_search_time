### Polymer Simulations with loop extruding factors and Experimental Comparison  
### *(For the paper: **Condensin Accelerates Long-Range Intra-Chromosomal Interactions**

This repository contains two components:

1. **`simulation/`** — 3D polymer MD simulation with Loop-Extruding Factors (LEFs)  
2. **`analysis/`** — A notebook + experimental data files used to convert simulation units to experimental units and compare simulation predictions of search time and distance between distal loci to experimental CICI data and average probe distance.

### simulation/
simulation.py : 
This script runs a 3D polymer molecular dynamics (MD) simulation with loop-extruding factors (LEFs) using polychrom and a LEFSimulator.
The 1D LEF simulation uses the lef-cython library (https://github.com/mirnylab/lefs-cython/tree/main).
The 3D simulation uses polychrom, a python wrapper of openmm, tailored for chromatin polymer simulations: https://github.com/open2c/polychrom

All parameters are specified in  parameters.py and the simulation uses helper functions from helper_functions.py.

Probes vs. all monomers
We distinguish between:

Probe monomers:

	•	A subset of monomers representing the LacO and TetO arrays at high temporal resolution
	•	Defined by: probe_indices = np.arange(0, N*M-1, 180)

All monomers:

	•	The full polymer configuration (all N*M beads), sampled at lower temporal resolution
	•	Used to compute Hi-C–like contact maps and other equilibrium observables.

Output files:
All output is stored in a folder named:
processivity_<processivity>_d_<d>/

### analysis/

The analysis/ folder contains a
- two .npy files: 1) experimental MSD for pair 4 +AID 2) corresponding lag times. 
- a notebook 

The analysis notebook takes precomputed probe-only trajectories from:

- `no_extrusion_precomputed/`  (no LEFs)
- `processivity_210_d_2100_precomputed/`  (with LEFs $\lambda=210kb$ , $d=2100kb$)

These folders must be downloaded from Zenodo (10.5281/zenodo.17777768) and placed inside the `analysis/` directory before running the notebook.

and performs three steps:

1. **Calibration to experimental units**  
2. **Inferring reaction radius**  
3. **Quantitative comparison: +AID vs −AID**