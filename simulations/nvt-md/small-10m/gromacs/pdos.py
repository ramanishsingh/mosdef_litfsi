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
import fileinput

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
        print(
            "The combined trajectory has {} frames = {} ps ".format(
                comb_traj.n_frames, comb_traj.n_frames * 5 / 1000
            )
        )
        comb_traj.save('comb_traj.h5')
        # print('All residues: %s' % [residue for residue in comb_traj.topology.residues])
        num_frames=500
        idxs=np.round(np.linspace(0, comb_traj.n_frames - 1, num_frames)).astype(int) 
        
        for idx in idxs:
            os.mkdir(str(idx))
            os.chdir(str(idx))
            shutil.copyfile("../../../pdos_files/submit_pdos.sh","submit_pdos.sh")
            shutil.copyfile("../../../pdos_files/litfsi_pdos.inp","litfsi_pdos.inp")
            shutil.copyfile("../../../pdos_files/pdos.py","pdos.py")
            shutil.copyfile("../../../pdos_files/get-smearing-pdos.py","get-smearing-pdos.py")
            shutil.copyfile("../../../pdos_files/post_process_pdos.py","post_process_pdos.py")
            shutil.copyfile("../../../pdos_files/energy_values_array.txt","energy_values_array.txt")
            shutil.copyfile("../../../pdos_files/new_pdos_1.py","new_pdos_1.py")
            frame=comb_traj[idx].save('litfsi.xyz')
            with fileinput.FileInput("litfsi_pdos.inp", inplace=True) as file:
                for line in file:
                    print(line.replace("LENGTH", str(length[0])), end='')
            with fileinput.FileInput("submit_pdos.sh", inplace=True) as file:
                for line in file:
                    print(line.replace("idx", str(idx)), end='')
            os.system("sbatch submit_pdos.sh")
            os.chdir("..")
        os.chdir("..")

 
if __name__ == "__main__":
    main()
