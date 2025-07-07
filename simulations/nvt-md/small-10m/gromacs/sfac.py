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
    data_path = "sfac_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    top = "../../10m_small.pdb"
    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T", "Functional")):
        print(temp, functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list = []
        for job in group:
            seed = job.sp.Seed
            length = job.sp.L
            trj_file = os.path.join(job.workspace(), "sample.xtc")
            full_traj = md.load(trj_file, top = top)
            if full_traj.n_frames>500:
                traj_list.append(full_traj) # discarding first 0 ps

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
            "The combined trajectory has {} frames".format(
                comb_traj.n_frames
            )
        )
        # print('All residues: %s' % [residue for residue in comb_traj.topology.residues])

        # First calculating the overall SFAC

        q, s = scattering.structure_factor(
            comb_traj, Q_range=(L_min, 200), n_points=2000, form="cromer-mann"
        )
        q = q / 10
        # s=s+1
        fig = plt.figure()
        plt.plot(q, s)
        q = np.reshape(q, (q.shape[0], 1))
        s = np.reshape(s, (q.shape[0], 1))
        c = np.concatenate((q, s+1), axis=1)
        np.savetxt("sfac.txt", c, header="q (1/AA)   S(q)")
        plt.grid(alpha=0.3)
        plt.xlabel("$\it{q} (1/\AA)$")
        plt.ylabel("$S(q)$")
        plt.savefig("sfac.pdf")
        plt.close()
        q_overall = q
        s_overall = s

        partials = scattering.structure_factor(
            comb_traj,
            Q_range=(L_min, 200),
            n_points=2000,
            form="cromer-mann",
            partial=True,
        )

        # print(partials.keys())

        # We have to calculate 7 sfacs, s_Li_Li, s_Li_TFSI, s_Li_water, s_TFSI_TFSI, s_TFSI_water, s_water_water, s_summed

        # Getting s_Li_Li
        s_Li_Li = partials[("Li1", "Li1")]
        plt.plot(q_overall, s_Li_Li, label="{}-{}".format("Li", "Li"))

        # Writing out the atom names in each of the residues
        TFSI_atoms = [
            "O1",
            "F1",
            "C1",
            "S1",
            "O2",
            "F2",
            "N1",
            "O3",
            "F3",
            "S2",
            "O4",
            "F4",
            "F5",
            "C2",
            "F6",
        ]
        water_atoms = ["O5", "H1", "H2"]
        Li_atoms = ["Li1"]

        # Making pair lists, and getting partial sfac

        # TFSI-TFSI

        TFSI_pairs = []

        for (atom1, atom2) in it.product(TFSI_atoms, repeat=2):
            TFSI_pairs.append((atom1, atom2))
        s_TFSI_TFSI = 0
        for pair in TFSI_pairs:
            s_TFSI_TFSI += partials[pair]

        # H2O-H20
        water_pairs = []

        for (atom1, atom2) in it.product(water_atoms, repeat=2):
            water_pairs.append((atom1, atom2))
        s_water_water = 0
        for pair in water_pairs:
            s_water_water += partials[pair]

        # Li-H2O
        Li_water_pairs = [("Li1", "O5"), ("Li1", "H1"), ("O5", "Li1"), ("H1", "Li1")]
        s_Li_water = 0
        for pair in Li_water_pairs:
            s_Li_water += partials[pair]

        # TFSI-H2O
        s_TFSI_water = 0
        TFSI_water_pairs = []
        for (atom1, atom2) in it.product(water_atoms, TFSI_atoms):
            TFSI_water_pairs.append((atom1, atom2))
        for (atom1, atom2) in it.product(TFSI_atoms, water_atoms):
            TFSI_water_pairs.append((atom1, atom2))
        for pair in TFSI_water_pairs:
            s_TFSI_water += partials[pair]

        # Li-TFSI
        s_Li_TFSI = 0
        Li_TFSI_pairs = []
        for (atom1, atom2) in it.product(Li_atoms, TFSI_atoms):
            Li_TFSI_pairs.append((atom1, atom2))
        for (atom1, atom2) in it.product(TFSI_atoms, Li_atoms):
            Li_TFSI_pairs.append((atom1, atom2))
        for pair in Li_TFSI_pairs:
            s_Li_TFSI += partials[pair]

        # Summed
        s_summed = 0

        for key in partials.keys():
            s_summed += partials[key]

        # plotting, already plotted Li-Li
        plt.grid(alpha=0.3)
        plt.xlim([0, 3])
        plt.xlabel("$\it{q} (1/\AA)$")
        plt.ylabel("$S(q)$")
        plt.plot(q_overall, s_summed, "o-", label="Summed")
        plt.plot(q_overall, s_overall, label="Total")
        plt.plot(q_overall, s_TFSI_TFSI, label="TFSI-TFSI")
        plt.plot(q_overall, s_water_water, label="H2O-H2O")
        plt.plot(q_overall, s_Li_water, label="Li-water")
        plt.plot(q_overall, s_TFSI_water, label="TFSI-water")
        plt.plot(q_overall, s_Li_TFSI, label="Li-TFSI")
        plt.legend()
        plt.savefig("decomposed.png")
        plt.close()
        s_dict = {}
        s_dict["Li-Li"] = s_Li_Li
        s_dict["TFSI-TFSI"] = s_TFSI_TFSI
        s_dict["H2O-H2O"] = s_water_water
        s_dict["Li-H2O"] = s_Li_water
        s_dict["TFSI-H2O"] = s_TFSI_water
        s_dict["Li-TFSI"] = s_Li_TFSI
        # Saving partial sfacs
        pairs = ["Li-Li", "TFSI-TFSI", "H2O-H2O", "Li-H2O", "TFSI-H2O", "Li-TFSI"]

        for pair in pairs:
            s = np.reshape(s_dict[pair], (q_overall.shape[0], 1))

            c = np.concatenate((q_overall, s), axis=1)
            np.savetxt("sfac_{}.txt".format(pair), c, header="q (1/AA)   S(q)")

        os.chdir("..")


if __name__ == "__main__":
    main()
