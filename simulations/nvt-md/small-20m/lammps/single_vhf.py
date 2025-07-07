import os
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import signac
import unyt as u
import shutil
from scipy.integrate import simps
import freud
import sys
from scattering.van_hove import compute_van_hove
from scattering.van_hove import compute_partial_van_hove
from scattering.utils.utils import get_dt
from scattering.utils.features import find_local_maxima
import argparse

def main(seed):

    # top file to be used for the project
    top="../../../20m_small.pdb"
    traj = md.load("mol_{}.lammpstrj".format(seed), top=top)
    print("traj has n_frames", traj.n_frames)
    traj = traj[10000:]
    chunk_length = 2002
    end_length = traj.n_frames
    start_frame=0
    vhf_mean=0
    windows_counted=0    
    skip_step = 1


    while start_frame<(end_length-chunk_length):
        windows_counted+=1
        end_frame=start_frame+chunk_length
        chunk = traj[start_frame:end_frame]
        print(f"Analyzing frames {start_frame} to {end_frame}...")
        chunk.time = np.linspace(0*1000/1000, (chunk_length-1)*1000/1000, len(chunk.time))
        dt = get_dt(chunk)
        print("dt is {}".format(dt))
        print('Starting vhf')
        r, g_r_t = compute_partial_van_hove(trj=chunk,
                                       chunk_length=chunk_length,
                                       selection1= "element Li",
                                       selection2= "element O and resname wat",
                                       r_range=(0, 0.75),
                                       n_bins=150,
                                       self_correlation=False)
        t = chunk.time

        vhf_mean=((windows_counted-1)*vhf_mean+g_r_t)/windows_counted
        start_frame+= skip_step
                


    vhf_mean = vhf_mean.T
    np.savetxt("vhf_{}.txt".format(seed), vhf_mean)
    print("The shape of one vhf is ", vhf_mean.shape)
    print("r", r)
    print("t", t)
 #       vhf_write  = vhf_mean[::1]
 #       t_write    = t[::1]


    np.savetxt(
        f"t.txt",
        t,
    )

    np.savetxt(
        f"r.txt",
        r,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run vhf for one seed')
    parser.add_argument("-seed", metavar='seed', required=True, help="seed")
    args = parser.parse_args()
    main(args.seed)
