import os
import freud
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import matplotlib
import shutil
import signac
import networkx as nx
from hbond_elliptical import calualateHBMap
import collections

def main():
    # All output data will be stored in rdf_data folder
    save_string_file = "" # This string will be saved to a file
    data_path = "hbond_elliptical_cluster_water"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    # top file to be used for the project
    top="../../20m_small.pdb"

    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T","Functional")):
        #Only computing RDFs for PBE functional
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

                traj_list.append(full_traj[5000:]) # discarding first 25 ps

        comb_traj=md.join(traj_list)
        # Add unit cell information
        comb_traj = md.Trajectory(
                comb_traj.xyz,
                comb_traj.top,
                unitcell_lengths = np.tile([length[0]/10, length[1]/10, length[2]/10], (comb_traj.n_frames,1)),
                unitcell_angles = np.tile([90.,90.,90.], (comb_traj.n_frames,1)),
            )
        print("The combined trajectory has {} frames = {} ps ".format(comb_traj.n_frames,comb_traj.n_frames*5/1000))
        save_string_file += "\n"+"The combined trajectory has {} frames = {} ps \n".format(comb_traj.n_frames,comb_traj.n_frames*5/1000)+"\n"
        box = freud.box.Box(Lx=length[0]/10,Ly= length[1]/10, Lz=length[2]/10)
        top = comb_traj.topology


        nbins_r = 400
        nbins_a = 400
        r_cutoff = 0.75
        skip_every_x_frames = 10
        sel_oxygen_head = 'name O5' ; sel_hydrogen = 'name H1 or name H2' ; sel_oxygen_tail = 'name O5'; list_names_hydrogen = ["H1", "H2"] ; list_names_oxygen_head = ["O5"] ; list_names_oxygen_tail = ["O5"] 
        rdf_output, inter_output, map_output,hbond,hbond_time = calualateHBMap(comb_traj, r_cutoff, nbins_r, nbins_a, skip_every_x_frames, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail)
     

        wat_O_indices=comb_traj.top.select("element O and resname wat")

        super_cluster_size = []
        total_water = len(wat_O_indices)

        for index in range(len(hbond)):
            print("==========")
            data = hbond[index]
            pairs = []
            for pair in data:
        
                pairs.append((pair[0], pair[2]))
            G = nx.Graph()
            G.add_edges_from(pairs)
            summation = 0 
            for connected_component in nx.connected_components(G):
                super_cluster_size.append(len(connected_component))
                summation+=len(connected_component)
                print(len(connected_component))
            if summation < total_water:
                super_cluster_size.extend([1]*(total_water-summation))
        print(super_cluster_size)

        counter=collections.Counter(super_cluster_size)
        frequency_dict=dict(counter)
        frequency= np.array(list(counter.values()))
        #print("The sum of all freq for one frame is {}".format(np.sum(frequency)/comb_traj.n_frames))
        frac_cluster_size=np.array(list(counter.keys()))
        p = frac_cluster_size.argsort()
        frac_cluster_size=frac_cluster_size[p]/total_water
        frequency=frequency[p]
        #print(frequency)
        prob=[]
        for i in range(frequency.shape[0]):
         #   print("i is {}".format(i))
         #   print("local frq is {}".format(frequency[i]))
         #   print("sum of all freq is {}".format(np.sum(frequency)))
         #   print("ratio is {}".format(frequency[i]/np.sum(frequency)))
            prob.append(frequency[i]*frac_cluster_size[i]/(len(hbond)))
        prob=np.array(prob)
        cum_prob=np.cumsum(prob)
        #print("prob is  {} ".format(prob))
        #print("cum prob is  {} ".format(cum_prob))
        #print("frac cluster size is {}".format(frac_cluster_size))
        plt.figure()
        plt.plot(frac_cluster_size,cum_prob)
        plt.grid(alpha=0.2)
        plt.xlim([0.0001, 1])
        plt.ylim([1e-6, 1])
        plt.xlabel("N_{water in cluster}/N_{total water in box}")
        plt.ylabel("Cumulative probability")
        plt.savefig("cum_prob_O_O.png")
        np.savetxt(
            f"cum_prob_O_O.txt",
            np.transpose(np.vstack([frac_cluster_size, cum_prob])),
            header="frac_cluster_size\tcum_prob",
        )

        plt.close()



        #####################################################################
        # Compute Ow–Ow neighbor distribution using hbond criterion
        #####################################################################
        neighbor_counts_all = []

        for hbonds_frame in hbond:  # Loop over each frame
            # Store which Ow each Ow H-bonds with in this frame
            neighbor_map = collections.defaultdict(set)

            for (o_tail, h, o_head) in hbonds_frame:
                # Skip if it's not a WAT oxygen
                #if top.atom(o_tail).resname != 'wat' or top.atom(o_head).resname != 'wat':
                #    continue
                # Add both directions to ensure symmetry (optional)
                neighbor_map[o_tail].add(o_head)
                neighbor_map[o_head].add(o_tail)

            # Count how many unique H-bonded Ow neighbors each Ow has
            for ow_idx in wat_O_indices:
                num_neighbors = len(neighbor_map[ow_idx])
                neighbor_counts_all.append(num_neighbors)

        # Make distribution
        unique_counts = np.unique(neighbor_counts_all)
        frequencies = np.array([neighbor_counts_all.count(k) for k in unique_counts])
        probabilities = frequencies / np.sum(frequencies)

        # Plotting
        plt.figure()
        plt.bar(unique_counts, probabilities, width=0.8, align='center')
        plt.xlabel("Number of H-bonded Ow neighbors")
        plt.ylabel("Probability")
        plt.grid(alpha=0.2)
        plt.savefig("Ow_Ow_hbond_neighbor_distribution.png")

        # Save data
        np.savetxt("Ow_Ow_hbond_neighbor_distribution.txt",
                   np.column_stack([unique_counts, probabilities]),
                   header="n_neighbors\tprobability")
        plt.close()


        os.chdir("..")

        text_file = open("hbond_analysis_info_{}_{}.txt".format(temp, functional), "w")
        n = text_file.write(save_string_file)
        text_file.close()

if __name__ == "__main__":
    main()

