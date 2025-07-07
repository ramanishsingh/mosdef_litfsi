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
from scattering.utils.features import minima_in_range
from scipy import optimize


def time_and_distance(chunk, unwrap_chunk, Li_indices, wat_O_indices, box, Li_O_water_minima):
    neighbors_Li_original = []
    neighbors_Li=[]
    num_water_neighbors_list = []
    for frame in range(chunk.n_frames):
        water_pairs = []
        points_Li=chunk.xyz[frame][Li_indices]
        points_O=chunk.xyz[frame][wat_O_indices]
        aq = freud.locality.AABBQuery(box, points_O)
        for point in range(points_Li.shape[0]):
        #for point in range(1):
            neighbor_list=[]
            query_result = aq.query(points_Li[point], dict(r_max=Li_O_water_minima/10))
            nlist = query_result.toNeighborList()
            num_water_neighbors_list.append(len(nlist))
            for pair in list(nlist):
                neighbor_list.append(pair[1])
            water_pairs.append(neighbor_list)
        neighbors_Li.append(water_pairs)
    neighbors_Li_original = neighbors_Li[0]
    print(neighbors_Li_original)
    print("Li has on avg this many wtaer neighbors", np.mean(num_water_neighbors_list))
    data = []
    for li_number, li_index in enumerate(Li_indices):
        if len(neighbors_Li_original[li_number])==0:
            print("0 initial neighbors")
            continue
        original_neighbors = neighbors_Li_original[li_number]
        saw=False
        for frame in range(chunk.n_frames):
            fraction = (len([i for i in neighbors_Li_original[li_number] if i in neighbors_Li[frame][li_number] ]))/(len(neighbors_Li_original[li_number]))
            if fraction<=0.501:
                distance = 10*np.linalg.norm(unwrap_chunk.xyz[frame][li_index]-unwrap_chunk.xyz[0][li_index])
                time = frame*1
                data.append([distance,time])
                saw = True
                break
        if not saw:
            print("did not see fraction 0.5")
    print(data)
    print(len(data))
    return data


def main():
    # All output data will be stored in vhf_data folder
    save_string_file = "" # This string will be saved to a file
    data_path = "li_o_neighbor_dist_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    # top file to be used for the project
    top="../../../20m_small.pdb"
    project = signac.get_project()
    for (temp), group in project.groupby(("T")): 
        functional = "MD"
        print(temp,functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list=[]
        unwrap_traj_list = []

        traj_lengths = []
        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            
            trj_file = os.path.join(job.workspace(), "sample.xtc")
            full_traj = md.load(trj_file, top=top)
            print("full traj has n_frames", full_traj.n_frames)
            traj_list.append(full_traj) 


            trj_file = os.path.join(job.workspace(), "sample_unwrapped.xtc")
            full_traj = md.load(trj_file, top=top)
            unwrap_traj_list.append(full_traj) 

            traj_lengths.append(full_traj.n_frames)


        print("All trajs loaded")
        print("There are {} trajs".format(len(traj_list)))
        box = freud.box.Box(Lx=length[0]/10,Ly= length[1]/10, Lz=length[2]/10)

        data_overall=[]
        for traj_number, traj in enumerate(traj_list):
            unwrap_traj = unwrap_traj_list[traj_number]
            r_min=0.01
            r_max=min(length)/20-0.01
            n_bins=150
            print("r_max is {} nm".format(r_max))

            pair_indices_1=traj.top.select("element Li")
            pair_indices_2=traj.top.select("resname wat and element O")
            gr = 0
            freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
            for a,b,c in zip(np.asarray(traj.unitcell_vectors), traj.xyz[:, pair_indices_1, :],traj.xyz[:, pair_indices_2, :]):
                freud_rdf.compute((a,b),query_points=c, reset=False)
            r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
            Li_O_water_minima=minima_in_range(r*10, gr, 2,4)[0]
            print("The first shell minima of Li-O(water) lies at {} \AA".format(Li_O_water_minima))




            wat_O_indices=traj.top.select("element O and resname wat")
            Li_indices=traj.top.select("element Li")
            TFSI_O_indices = traj.top.select("resname TF2 and element O")
            TFSI_S_indices = traj.top.select("resname TF2 and element S")
            chunk_length = min(traj_lengths)-39000
            print("chunk length being used is ", chunk_length, traj_lengths)
            end_length = traj.n_frames
            start_frame=0
            windows_counted=0
            skip_step = 250

            while start_frame<(end_length-chunk_length):
                windows_counted+=1
                end_frame=start_frame+chunk_length
                chunk = traj[start_frame:end_frame]
                unwrap_chunk = unwrap_traj[start_frame:end_frame]
                print(f"Analyzing frames {start_frame} to {end_frame}...")
                chunk.time = np.linspace(0*1000/1000, (chunk_length-1)*5/1000, len(chunk.time))
                dt = get_dt(chunk)
                data = time_and_distance(chunk, unwrap_chunk, Li_indices, wat_O_indices, box, Li_O_water_minima)
                print(data)
                data_overall.extend(data)
                start_frame+= skip_step
        # Saving all data to a text file
        with open("all_data.txt", "w") as file:
            for item in data_overall:
                file.write("{} {}\n".format(item[0], item[1]))

        # Plotting 2D heatmap
        data_overall = np.array(data_overall)
        distances = data_overall[:, 0]
        times = data_overall[:, 1]
        plt.figure(figsize=(8, 6))
        plt.hist2d(distances, times, bins=(50, 50), cmap=plt.cm.jet)
        plt.colorbar(label='Counts')
        plt.xlabel('Distance (Angstroms)')
        plt.ylabel('Time (ps)')
        plt.title('2D Heatmap of Distance vs. Time')
        plt.savefig("heatmap_dist_vs_time.pdf")
        plt.close()

        os.chdir("..")


if __name__ == "__main__":
    main()

