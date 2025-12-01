import numpy as np

Nsteps = 1_000_000_000 # total simulation steps, all simulation in figures ran for 1_000_000_000 steps
Nthermalization = 5_000_000
period = 100 # save every period steps for the FPT calculation (fast acquisition rate)
periodall = 100000 # save every periodall steps for contact frequency calculation (slower acquisition rate)
restartSimulationEveryLEFstep = 20

N = 1800
M = 10
probe_indices = np.arange(0, N*M-1, 180) # write check for probe_indices vs N

col = 0.2 
steps = 10000 # steps between two LEF updates. For steps=5000 sets extrusion speed to 1 kb/s for symetric LEFs
GPU_ID = '0'

processivity = 210
d = 2100

confinement_radius = 50

# parameters for smc bonds
smcBondWiggleDist = 0.2
smcBondDist = 0.5
