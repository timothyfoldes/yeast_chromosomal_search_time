LEF–Polymer Simulation with Two Acquisition Rates

This script runs a 3D polymer molecular dynamics (MD) simulation with loop-extruding factors (LEFs) using polychrom and a custom LEFSimulator. It is designed to produce trajectories suitable for:
	1.	Fast-acquisition probe trajectories for first-passage time (FPT) / MS2-like analyses
	2.	Slow-acquisition full-polymer trajectories for equilibrium properties, such as Hi-C–like contact maps

All parameters are specified at the top of the script (or in parameters.py), and the simulation runs on a single CUDA GPU.

⸻

Model overview
	•	Polymer:
	•	Total monomers: N * M
	•	M independent polymers (chromosomes), each with N monomers
	•	Here: N = 1800, M = 10 (10 polymers of 1800 beads each)
	•	Geometry & forces:
	•	Polymers confined in a 3D sphere of radius confinement_radius
	•	Harmonic bonds along the backbone, no bending rigidity (k = 0 for angles)
	•	Short-range repulsive non-bonded potential between all non-bonded pairs
	•	Loop extrusion:
	•	LEFs are simulated in 1D along the polymer using LEFSimulator
	•	A total of N_LEFS = (M * N) / d extruders, with processivity processivity
	•	LEFs load everywhere along the polymer, unload with a rate set by processivity,
and are forced to release at telomere indices so that extruders cannot cross between chromosomes
	•	The resulting LEF positions at each LEF step are stored and then used to dynamically update SMC-like bonds in the 3D MD simulation via bondUpdater

⸻

Time scales and acquisition

There are three key time/step scales:
	•	Nsteps – total number of MD steps per run
	•	steps – number of MD steps between LEF updates
	•	LEF_steps = Nsteps // steps – total number of LEF configurations
	•	restartSimulationEveryLEFstep – how many LEF updates are grouped into one MD run before restarting the simulation object
	•	period – fast sampling interval (in MD steps)
	•	Used for FPT calculations, where we want high temporal resolution for a few loci
	•	periodall – slow sampling interval (in MD steps)
	•	Used for equilibrium properties (e.g. contact maps), where we can sample less often but need all monomers

The checks at the top (asserts) ensure that all these quantities are compatible (e.g. steps is a multiple of period, etc.).

⸻

Probes vs. all monomers

We distinguish between:
	•	Probe monomers:
	•	A subset of monomers used as “reporter loci”, sampled at high temporal resolution
	•	Defined by:

probe_indices = np.arange(0, N*M-1, 180)

i.e. every 180th monomer along the concatenated polymers.

	•	This spacing corresponds to the genomic distance between the LacO and TetO arrays (“probes”) in the experimental system.

	•	All monomers:
	•	The full polymer configuration (all N*M beads), sampled at lower temporal resolution
	•	Used to compute Hi-C–like contact maps and other equilibrium observables.

⸻

Simulation workflow
	1.	Parameter setup
	•	Reads Nsteps, Nthermalization, N, M, steps, period, periodall, restartSimulationEveryLEFstep, processivity, d, col, confinement_radius, SMC bond parameters, etc.
	•	Computes:

LEF_steps = Nsteps // steps
N_LEFS = (M * N) // d


	•	Sets up LEF loading/unloading probability arrays, including forced release at telomeric ends so extruders don’t cross chromosomes.

	2.	LEF simulation (1D)
	•	Initializes a LEFSimulator with N_LEFS extruders on a 1D lattice of length M * N.
	•	Runs a short thermalization (LEFs.steps(0, processivity * 10)) to equilibrate LEF distribution.
	•	Then runs LEF_steps iterations:

for k in range(LEF_steps):
    LEFs.steps(k, k + 1)
    LEFpositions[k, :, :] = LEFs.get_LEFs()


	•	Saves all LEF positions:

np.save(f"{folder}/LEFpositions.npy", LEFpositions)


	3.	Bond updater (“milker”)
	•	bondUpdater(LEFpositions) translates 1D LEF positions into 3D SMC-like bonds between monomers.
	•	At run time, it activates/deactivates these bonds in the MD simulation according to the precomputed LEF trajectories.
	4.	MD simulation (3D)
The simulation is run in chunks, restarting the Simulation object every restartSimulationEveryLEFstep LEF updates:
	•	Create Simulation (CUDA, variable Langevin integrator, collision rate col)
	•	Set initial positions using a cubic “grown” conformation:

data = grow_cubic(N*M, boxSize=int((N*M)**(1/3)*1.2))
sim.set_data(data)


	•	Add forces:
	•	Spherical confinement
	•	Polymer chains with harmonic bonds, angle forces, and repulsive non-bonded interactions
	•	Initialize SMC bonds via milker.setParams(...) and milker.setup(...)
	•	Energy minimization on the first iteration
	•	Thermalization on the first iteration:
	•	High-friction run for Nthermalization steps
	•	Short relaxation at nominal collision rate
	•	Production run:
	•	For each LEF block, run restartSimulationEveryLEFstep * (steps // period) small MD blocks of length period:
	•	Record probe positions every period
	•	Record all positions every periodall
	•	Advance LEF bonds via milker.step(...) every steps (in MD steps)

⸻

Output files

All output is stored in a folder named:

processivity_<processivity>_d_<d>/

Within this folder:

1. LEF positions (1D)
	•	LEFpositions.npy
Shape: (LEF_steps, N_LEFS, 2)
	•	For each LEF step and each extruder, stores the left and right anchor positions along the 1D polymer.
	•	Used by bondUpdater to impose active SMC bonds in 3D.

2. Fast-acquisition probe trajectories (for FPT calculations)

For each simulation chunk iteration:
	•	trajprobe_<iteration>.npy
	•	3D positions of probe monomers only at fast acquisition rate period.
	•	Shape: (N_frames_probe, N_probes, 3)
	•	N_probes = len(probe_indices) (here one probe every 180 monomers)
	•	These are intended for computing:
	•	First-passage times (FPT) to a given distance threshold between specific probes
	•	MS2-like transcriptional readouts based on distance-dependent contact events
	•	timeprobe_<iteration>.npy
	•	Simulation times (same length as trajprobe) in simulation time units.
	•	Already corrected for the cumulative time across iterations.

Use case:
Compute distances between specific probe pairs (e.g. those mimicking LacO/TetO arrays) and quantify:
	•	Mean FPT
	


3. Slow-acquisition full trajectories (for equilibrium / Hi-C–like statistics)

For each simulation chunk iteration:
	•	trajall_<iteration>.npy
	•	3D positions of all monomers sampled every periodall MD steps.
	•	Shape: (N_frames_all, N * M, 3)
	•	Intended for equilibrium observables:
	•	Hi-C–like contact maps
	•	Distance vs genomic separation (R(s))
	•	Average spatial distances and other static statistics
	•	timeall_<iteration>.npy
	•	Simulation times corresponding to the frames in trajall_<iteration>.npy.

4. Reporter & completion flag
	•	HDF5 files written by HDF5Reporter (optional, mainly for visualization/debugging with polychrom tools).
	•	completed.txt
	•	A simple flag file containing "completed", written at the end of the simulation to indicate successful completion.

⸻

Typical downstream analysis
	•	First-passage times (FPT)
Use trajprobe_* and timeprobe_*:
	1.	Select the two probe indices corresponding to the loci of interest (e.g. LacO and TetO analogues).
	2.	Compute their 3D distance vs time.
	3.	From this 1D distance signal, compute FPTs to a chosen distance cutoff (contact threshold), possibly averaging over multiple probe pairs and chunks.
	•	Hi-C–like contact maps and equilibrium properties
Use trajall_* and timeall_*:
	1.	Bin monomers along the polymer (if desired).
	2.	For each saved frame, mark pairs of monomers that are closer than a chosen cutoff.
	3.	Average over frames to obtain a contact probability matrix.
	4.	Compute derived observables like P(s), R(s), etc.

⸻

This setup thus cleanly separates:
	•	Fast, sparse sampling of probes → suited for kinetic quantities (FPT, on/off times)
	•	Slow, dense sampling of all monomers → suited for equilibrium quantities (Hi-C, R(s), static structure)

while both are generated from the same underlying LEF–polymer dynamics.