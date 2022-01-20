import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import signac
import shutil
from scipy import stats
import freud
import mbuild as mb
import collections


def unique(list1):

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list


def main():
    data_path = "cumulative_prob_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)
    
    top="../../20m_small.pdb"
    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T","Functional")):
        if functional =="BLYP":
            continue
        print(temp,functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list=[]
        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            full_traj = md.load(job.fn("litfsi-pos-1.xyz"), top=top)
            if full_traj.n_frames>5000:

                traj_list.append(full_traj[5000:])
        comb_traj=md.join(traj_list)
        # Add unit cell information
        comb_traj = md.Trajectory(
                comb_traj.xyz,
                comb_traj.top,
                unitcell_lengths = np.tile([length[0]/10, length[1]/10, length[2]/10], (comb_traj.n_frames,1)),
                unitcell_angles = np.tile([90.,90.,90.], (comb_traj.n_frames,1)),
            )
        #print(comb_traj)
        print("The combined trajectory has {} frames = {} ps ".format(comb_traj.n_frames,comb_traj.n_frames*5/1000))
        #print('All residues: %s' % [residue for residue in comb_traj.topology.residues])
        box = freud.box.Box(Lx=length[0]/10,Ly= length[1]/10, Lz=length[2]/10)
        num_residues=comb_traj.n_residues

        # com_traj  is the COM traj and comb_traj is the complete traj, now find the rdfs
        r_min=0.01
        r_max=min(length)/20-0.01
        n_bins=150
        print("r_max is {} nm".format(r_max))
       


        wat_O_indices=comb_traj.top.select("element O and resname wat")

        super_cluster_size=[]
        super_num_clusters=[]
        for frame in range(comb_traj.n_frames):
            num_neighbors=[]
            points_O=comb_traj.xyz[frame][wat_O_indices]
            points_O=box.wrap(points_O)
            system = freud.AABBQuery(box, points_O)
            cl = freud.cluster.Cluster()
            cl.compute(system, neighbors={"r_max": 0.33})
            for cluster_id in range(cl.num_clusters):
                cluster_system = freud.AABBQuery(
                                                 system.box, system.points[cl.cluster_keys[cluster_id]]
                                                 )
                super_cluster_size.append(len(cl.cluster_keys[cluster_id]))
        counter=collections.Counter(super_cluster_size)
        frequency_dict=dict(counter)
        frequency= np.array(list(counter.values()))
        print("The sum of all freq for one frame is {}".format(np.sum(frequency)/comb_traj.n_frames))
        frac_cluster_size=np.array(list(counter.keys()))
        p = frac_cluster_size.argsort() 
        frac_cluster_size=frac_cluster_size[p]/points_O.shape[0]
        frequency=frequency[p]
        print(frequency)
        prob=[]
        for i in range(frequency.shape[0]):
            print("i is {}".format(i))
            print("local frq is {}".format(frequency[i]))
            print("sum of all freq is {}".format(np.sum(frequency)))
            print("ratio is {}".format(frequency[i]/np.sum(frequency)))
            prob.append(frequency[i]/np.sum(frequency))
        prob=np.array(prob)
        cum_prob=np.cumsum(prob)
        print("prob is  {} ".format(prob))
        print("cum prob is  {} ".format(cum_prob))
        print("frac cluster size is {}".format(frac_cluster_size))
        plt.figure()
        plt.loglog(frac_cluster_size,cum_prob)
        plt.grid(alpha=0.2)
        plt.xlim([0.0001, 1])
        plt.ylim([1e-6, 1])
        plt.savefig("cum_prob_O_O.png")
        np.savetxt(
            f"cum_prob_O_O.txt",
            np.transpose(np.vstack([frac_cluster_size, cum_prob])),
            header="frac_cluster_size\tcum_prob",
        )

        plt.close()


        TFSI_S_indices=comb_traj.top.select("element S and resname TF2")

        super_cluster_size=[]
        super_num_clusters=[]
        for frame in range(comb_traj.n_frames):
            num_neighbors=[]
            points_S=comb_traj.xyz[frame][TFSI_S_indices]
            points_S=box.wrap(points_S)
            system = freud.AABBQuery(box, points_S)
            cl = freud.cluster.Cluster()
            cl.compute(system, neighbors={"r_max": 0.56})
            for cluster_id in range(cl.num_clusters):
                cluster_system = freud.AABBQuery(
                                                 system.box, system.points[cl.cluster_keys[cluster_id]]
                                                 )
                super_cluster_size.append(len(cl.cluster_keys[cluster_id]))
        counter=collections.Counter(super_cluster_size)
        frequency_dict=dict(counter)
        frequency= np.array(list(counter.values()))
        print("The sum of all freq for one frame is {}".format(np.sum(frequency)/comb_traj.n_frames))
        frac_cluster_size=np.array(list(counter.keys()))
        p = frac_cluster_size.argsort()
        frac_cluster_size=frac_cluster_size[p]/points_S.shape[0]
        frequency=frequency[p]
        print(frequency)
        prob=[]
        for i in range(frequency.shape[0]):
            prob.append(frequency[i]/np.sum(frequency))
        prob=np.array(prob)
        cum_prob=np.cumsum(prob)
        print("prob is  {} ".format(prob))
        print("cum prob is  {} ".format(cum_prob))
        print("frac cluster size is {}".format(frac_cluster_size))
        plt.figure()
        plt.loglog(frac_cluster_size,cum_prob)
        plt.grid(alpha=0.2)
        plt.xlim([0.0001, 1])
        plt.ylim([1e-6, 1])
        plt.savefig("cum_prob_S_S.png")
        np.savetxt(
            f"cum_prob_S_S.txt",
            np.transpose(np.vstack([frac_cluster_size, cum_prob])),
            header="frac_cluster_size\tcum_prob",
        )

        plt.close()


 
        os.chdir("..")    
if __name__ == "__main__":
    main()
