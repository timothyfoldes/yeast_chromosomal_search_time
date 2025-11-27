import os
import sys
import time
import numpy as np
from polychrom import forcekits, forces, simulation
from polychrom.hdf5_format import HDF5Reporter
from polykit.generators.initial_conformations import grow_cubic
from lefs_cython.simple import LEFSimulator
from helper_functions import bondUpdater
from parameters import *

LEF_steps = Nsteps//steps

sys.stderr.write(f'Running on GPU {GPU_ID} \n')
sys.stderr.write(f"Will do {Nsteps} simulation steps \n")
sys.stderr.write(f"will reconstruct every {(LEF_steps//restartSimulationEveryLEFstep)*steps} \n") #simInitsTotal =  (LEF_steps) // restartSimulationEveryLEFstep = (Nsteps//steps) // restartSimulationEveryLEFstep

### Setting up probability arrays for LEF simulation

N_LEFS = (M*N)//d

load_array = 1 * np.ones(((M*N)))
unload_array = np.ones(((M*N), 5))/(0.5*processivity)  
capture_array = np.zeros(((M*N), 2))  # no CTCF
release_array = np.zeros((M*N))
pause_array = np.zeros((M*N))  # no pausing

telomer_indices = [(i+1) * N-1 for i in range(M)] + [i * N for i in range(M)]
unload_array[telomer_indices] = 2 # release at ends of polymer so that LEFs don't cross chromosomes (=2 because release probability is averaged over both legs)

folder = 'processivity_' + str(processivity) + '_d_' + str(d)
sys.stderr.write(f"running processivity = {processivity}, d= {d}\n")

### GENERATE LEF POSITIONS

LEFpositions = np.zeros((LEF_steps, N_LEFS, 2))
LEFs = LEFSimulator(N_LEFS, (M*N), load_array, unload_array, capture_array, release_array, pause_array, skip_load=False)
LEFs.steps(0,processivity*10) # thermalize LEFs
for k in range(LEF_steps):
    LEFs.steps(k,k+1)
    LEFpositions[k,:,:] = LEFs.get_LEFs()

LEFpositions = np.array(LEFpositions)

os.makedirs(folder, exist_ok=True)
np.save(f"{folder}/LEFpositions.npy", LEFpositions)

sys.stderr.write("Making milker \n")
milker = bondUpdater(LEFpositions)

sys.stderr.write("Making reporter \n")

reporter = HDF5Reporter(folder=folder, max_data_length=100, overwrite=True, blocks_only=False)

### INITIAL POSITIONS

data = grow_cubic(N*M,boxSize=int((N*M)**(1/3)*1.2))
sys.stderr.write("EXTRUSION IS ON \n")

### CHECKS

assert ((restartSimulationEveryLEFstep * steps) % period) == 0
assert (periodall % period) == 0
assert (steps % period) == 0
assert (Nsteps % period) == 0
assert (periodall % steps) == 0

counter = 0

for iteration in range((Nsteps//steps) // restartSimulationEveryLEFstep): # Nsteps // steps = LEFsteps
    sys.stderr.write("Starting simulation \n")

    ### SIMULATION OBJECT CREATION

    sim = simulation.Simulation(
        platform="CUDA",
        integrator="variableLangevin",
        error_tol=0.003,
        GPU=GPU_ID,
        collision_rate=col,
        N=N*M,
        save_decimals=2,
        PBCbox=False,
        reporters=[reporter],
        verbose=False
    )  
    sys.stderr.write("Simulation object created \n")
    ### INITIAL POSITIONS
    
    sim.set_data(data) # set initial positions

    ### FORCES
    # Confinement
    sim.add_force(forces.spherical_confinement(sim,r=confinement_radius, k=1))
    sys.stderr.write(f"adding confinement with radius {confinement_radius}\n")
    
    # polymer backbone
    sim.add_force(
        forcekits.polymer_chains(
            sim,
            chains=[(i * N, (i + 1) * N, False) for i in range(M)],
            bond_force_func=forces.harmonic_bonds,
            bond_force_kwargs={
                "bondLength": 1,
                "bondWiggleDistance": 0.05,  # Bond distance will fluctuate +- 0.05 on average
            },
            angle_force_func=forces.angle_force,
            angle_force_kwargs={
                "k": 0., # no bending rigidity
            },
            
            nonbonded_force_func=forces.polynomial_repulsive,
            nonbonded_force_kwargs={
                "trunc": 3.0, 
            },
            except_bonds=True,
        )
    )

    ### MILKER 

    # ------------ initializing milker; adding bonds ---------
    # copied from addBond
    kbond = sim.kbondScalingFactor / (smcBondWiggleDist ** 2)
    bondDist = smcBondDist * sim.length_scale

    activeParams = {"length":bondDist,"k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    milker.setParams(activeParams, inactiveParams)
    
    # this step actually puts all bonds in and sets first bonds to be what they should be
    milker.setup(bondForce=sim.force_dict['harmonic_bonds'],
                blocks=restartSimulationEveryLEFstep)

    ### ENERGY MINIMIZATION

    if iteration==0:
        sim.local_energy_minimization() 
        previous_time_all = 0
        previous_time_probe = 0
    sim._apply_forces()

    ### THERMALIZATION

    if iteration==0:
        sys.stderr.write(f"Runing thermalization for {Nthermalization} steps \n")
        sim.integrator.setFriction = 0.03 # set collision rate to 0.01 for faster thermalization
        for i in range(Nthermalization//10_000):
            sim.do_block(10_000, save=False)
        
        sim.integrator.setFriction = col # set collision rate to back to normal
        sim.do_block(20_000, save=False) # run a bit more to let the velocities relax
        previous_time_all = 0
        previous_time_probe = 0
        therm_time = sim.time
    else:
        therm_time = 0


    ### SIMULATION

    trajall, trajprobe =  [], [] # to store trajectories
    timeall, timeprobe = [], []


    for i in range(restartSimulationEveryLEFstep*(steps//period)):
        
        ### RUN SIMULATION  
              
        if i%(10_000//period) == 0: # every 10_000 steps print sim state
            sim.do_block(steps=period, save=False, verbose=True)
            counter += 1
        else:
            sim.do_block(steps=period, save=False, verbose=False)
            counter += 1
        
        
        ### UPDATE BONDS

        posall = sim.get_data() # get all positions
        posprobes = np.array(posall[probe_indices]) # get only probe positions
        simtime = sim.time - therm_time #get simulation time in picoseconds

        ### GET PROBE POSITIONS

        if counter%1 == 0: # record probe positions every period
            trajprobe.append(posprobes)
            timeprobe.append(simtime)

        ### GET ALL POSITIONS

        if (counter % (periodall//period)) == 0: # record all positions every periodall
            trajall.append(posall)
            timeall.append(simtime)

        if (i < restartSimulationEveryLEFstep*(steps//period) - 1) and (counter % (steps//period)) == 0:
            curBonds, pastBonds = milker.step(sim.context)  # this updates bonds. You can do something with bonds here
    
    ### SAVE TRAJECTORIES

    timeall = np.array(timeall) + previous_time_all
    timeprobe = np.array(timeprobe) + previous_time_probe

    np.save(f"{folder}/trajall_{iteration}.npy", np.array(trajall))
    np.save(f"{folder}/trajprobe_{iteration}.npy", np.array(trajprobe))
    np.save(f"{folder}/timeall_{iteration}.npy", timeall)
    np.save(f"{folder}/timeprobe_{iteration}.npy", timeprobe)

    previous_time_all = timeall[-1]
    previous_time_probe = timeprobe[-1]

    ### 
    data = sim.get_data()  # save data and step, and delete the simulation
    del sim
    reporter.blocks_only = True  # Write output hdf5-files only for blocks
    time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)
    
reporter.dump_data()
with open(f"{folder}/completed.txt", "w") as f:
    f.write("completed")
time.sleep(1)