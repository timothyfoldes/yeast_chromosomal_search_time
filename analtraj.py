import freud
import numpy as np
from scipy import fft
from scipy.special import erfc
from scipy.optimize import curve_fit


def M2(t, Gamma2, tau, J):
    tau = (J/Gamma2)**2
    return 2 * Gamma2 * np.sqrt(t) * (1 - np.exp(-tau / (np.pi * t))) + 2 * J * erfc(np.sqrt(tau / (np.pi * t)))

def compute_msd(trajectories):
    """
    compute the mean squared displacement of a trajectory.
    Trajectores.shape must be (Nframes, Nparticles, 3)
    """

    box = freud.box.Box(Lx=10000, Ly=10000, Lz=10000)
    msd = freud.msd.MSD(box=box, mode='window')
    msd.compute(trajectories)
    return msd.msd

def compute_fpts(positions, target_size, Nsamples):
    # Create a boolean array where the condition is met
    target_reached = positions < target_size
        
    # Precompute first passage times for all positions
    first_passage_times = np.where(target_reached)[0]

    # Randomly choose Nsamples starting points after 100 to avoid the transient
    ts = np.random.randint(100, len(positions), size=Nsamples)

    # For each t, find the index in first_passage_times where first_passage_times >= t
    idxs = np.searchsorted(first_passage_times, ts, side='left')

    # Keep only indices where the target is reached after t
    valid = idxs < len(first_passage_times)
    ts = ts[valid]
    idxs = idxs[valid]

    # Calculate the first passage times
    fpts = first_passage_times[idxs] - ts

    if len(fpts) == 0:
        mean_fpts = len(positions)  # If target is never reached
    else:
        mean_fpts = np.mean(fpts)
    return np.array(fpts)


def extract_rij(traj, k, shift):
    """
    takes in trajectory AFTER being sliced into polymer chains, i.e (M, T, Nprobe, 3)
    then extracts rij for pairs of particles separated by k monomers
    output has dimension (number of pairs = M* (Nrobe//shift) - k, number of frames, 3). It is ready to be fed to compute_msd
    """
    Nprobe = traj.shape[-2]
    rij = np.array([(traj[:,:,i*shift,:] - traj[:,:,i*shift+k,:]) for i in range(0, (Nprobe-k)//shift)])
    return np.swapaxes(rij.reshape(-1, rij.shape[2], rij.shape[3]), 0, 1)

def compute_msd_two(traj_reshaped, k, shift):
    """
    computes the mean squared displacement of the trajectory
    """

    rij = extract_rij(traj_reshaped,k, shift)
    MSD_two = compute_msd(rij)
    rij2mean = np.mean(np.linalg.norm(rij, axis = -1)**2, axis = 0).mean(axis = 0)
    rij2mean = np.array(rij2mean)
    MSD_two = np.array(MSD_two)
    return MSD_two, rij2mean

def dct2(frames):
    N = np.shape(frames)[-2]
    inv_norm = 1 / (2 * N)
    return inv_norm * fft.dct(frames, axis=-2)

def modes2(frames):
    """
    frames.shape = (1024, 100, 3)
    where
    1024: frames number   (axis=-3 ou 0)
     100: signal size     (axis=-2 ou 1)
       3: space dimension (axis=-1 ou 2)
    """

    return np.square(np.abs(dct2(frames))).mean(axis=-3).sum(axis=-1)

def func(x, a, b):
    return a*x + b

def linfit(x, y, p, q):
    """
    fits a line to the data
    """
    param, pcov = curve_fit(func, x[p:q], y[p:q])
    return param

### DEFUNCT 

# def compute_fpts_old(positions, target_size, Nsamples = 10_000):
#     fpts = []
#     for _ in range(Nsamples):
#         t = np.random.randint(100, len(positions)) # start at 100 to avoid the initial transient
#         positions_conditionned = positions[t:]
#         if np.argmax(positions_conditionned < target_size).sum() == 0: # if the target is never reached pass
#             pass
#         else: # if the target is reached, record the first passage time
#             fpts.append(np.argmax(positions_conditionned < target_size))
#     if len(fpts) == 0:
#         mean_fpts =  len(positions)
#     else:
#         mean_fpts = np.mean(fpts)
#     return np.array(fpts)

def compute_fpts_old(positions, target_size, Nsamples):
    # Create a boolean array where the condition is met
    target_reached = positions < target_size
    
    # Precompute first passage times for all starting positions
    first_passage_times = np.where(target_reached)[0]

    fpts = []
    
    for _ in range(Nsamples):
        # Randomly choose a starting point after 100 to avoid the transient
        t = np.random.randint(100, len(positions))
        
        # Check if the target is ever reached after t
        if first_passage_times[first_passage_times >= t].size > 0:
            # Find the first index after t where the target is reached
            fpt = first_passage_times[first_passage_times >= t][0] - t
            fpts.append(fpt)
    
    if len(fpts) == 0:
        mean_fpts = len(positions)  # If target is never reached
    else:
        mean_fpts = np.mean(fpts)

    return np.array(fpts)
