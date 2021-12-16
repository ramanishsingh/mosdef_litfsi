import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import signac
import shutil
from scipy import stats
import freud
import mbuild as mb
from scattering import scattering
import itertools as it


def main():
    data_path = "pdos_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    top = "../../10m_small.pdb"
    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T", "Functional")):
        print(temp, functional)
        if functional=="BLYP":
            continue
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list = []
        for job in group:
            seed = job.sp.Seed
            length = job.sp.L
            full_traj = md.load(job.fn("litfsi-pos-1.xyz"), top=top)
            print(full_traj.n_frames)
            if full_traj.n_frames > 5000:

                traj_list.append(full_traj[5000:])

        comb_traj = md.join(traj_list)

        # Add unit cell information
        comb_traj = md.Trajectory(
            comb_traj.xyz,
            comb_traj.top,
            unitcell_lengths=np.tile(
                [length[0] / 10, length[1] / 10, length[2] / 10],
                (comb_traj.n_frames, 1),
            ),
            unitcell_angles=np.tile([90.0, 90.0, 90.0], (comb_traj.n_frames, 1)),
        )
        print(comb_traj)
        L_min = 2 * 10 * np.pi / length[0]
        print("L_min is {}".format(L_min))
        print(
            "The combined trajectory has {} frames = {} ps ".format(
                comb_traj.n_frames, comb_traj.n_frames * 5 / 1000
            )
       )
        comb_traj.save('comb_traj.h5')
        print("loading")
        comb_traj_new=md.load_hdf5('comb_traj.h5')
        print("loaded")

        os.chdir("..")

 
if __name__ == "__main__":
    main()
