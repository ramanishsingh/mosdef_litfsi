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

def fractions_time(chunk, Li_indices, wat_O_indices, box, Li_O_water_minima):
    neighbors_Li_original = []
    neighbors_Li=[]
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

            for pair in list(nlist):
                neighbor_list.append(pair[1])
            water_pairs.append(neighbor_list)
        neighbors_Li.append(water_pairs)
    neighbors_Li_original = neighbors_Li[0]

    #print("Li indices", Li_indices)
    #print("len neighbors li original", len(neighbors_Li_original))

    non_zero_li_indices = []
    for li_index in Li_indices:
        if len(neighbors_Li_original[li_index])>0:
            non_zero_li_indices.append(li_index)
    
    #print("neighbors_Li_original", neighbors_Li_original)
    #print("non_zero_li_indices", non_zero_li_indices)

    fraction = []
    for frame in range(len(neighbors_Li)):
        fraction_individual_list = [] 
        li_number = 0
        for li_index in non_zero_li_indices:
            li_number+=1
     #       print("This is li_number", li_number)
     #       print("and has neighbors =", len(neighbors_Li_original[li_index]))
            if len(neighbors_Li_original[li_index]) ==0:
                print("no neighbors found")
                return "No neighbors"
            fraction_individual =  (len([i for i in neighbors_Li_original[li_index] if i in neighbors_Li[frame][li_index] ]))/(len(neighbors_Li_original[li_index]))
            fraction_individual_list.append(fraction_individual)
        #print("fraction individual list", fraction_individual_list)
        fraction.append(np.mean(fraction_individual_list))
    return np.array(fraction)

def single_exp(t, tau):
    return np.exp(-t / tau)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = 0
    for idx in range(len(array)):
        if array[idx]<value:
            return idx, array[idx]

    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def main():
    # All output data will be stored in vhf_data folder
    save_string_file = "" # This string will be saved to a file
    data_path = "li_o_neighbor_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    # top file to be used for the project
    top="../../20m_small.pdb"
    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T", "Functional")):
        #Only computing VHFs for PBE functional
        if functional == "BLYP":
            continue
        print(temp,functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list=[]


        traj_lengths = []
        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            trj_file = os.path.join(job.workspace(), "litfsi-pos-1.xyz")
            full_traj = md.load(trj_file, top=top)
            print("full traj has n_frames", full_traj.n_frames)
            #if full_traj.n_frames<25000:
            #    continue
            if full_traj.n_frames>5000:
                traj = md.Trajectory(
                        full_traj.xyz,
                        full_traj.top,
                        unitcell_lengths = np.tile([length[0]/10, length[1]/10, length[2]/10], (full_traj.n_frames,1)),
                        unitcell_angles = np.tile([90.,90.,90.], (full_traj.n_frames,1)),
                    )
                traj_list.append(traj[5000:])
                traj_lengths.append(full_traj.n_frames-5000)

        #get the shortest traj and reduce all trajs to that length
        print("All trajs loaded")
        print("There are {} trajs".format(len(traj_list)))
        #print(traj_lengths);os.chdir("..");continue
        box = freud.box.Box(Lx=length[0]/10,Ly= length[1]/10, Lz=length[2]/10)


        fraction_overall_list = []
        for traj in traj_list:
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
            chunk_length = min(traj_lengths)-2000
            print("chunk length being used is ", chunk_length, traj_lengths)
            end_length = traj.n_frames
            start_frame=0
            windows_counted=0
            skip_step = 500
            fractions= []
            fractions_mean = 0

            while start_frame<(end_length-chunk_length):
                windows_counted+=1
                end_frame=start_frame+chunk_length
                chunk = traj[start_frame:end_frame]
                print(f"Analyzing frames {start_frame} to {end_frame}...")
                chunk.time = np.linspace(0*1000/1000, (chunk_length-1)*5/1000, len(chunk.time))
                dt = get_dt(chunk)
                fraction = fractions_time(chunk, Li_indices, wat_O_indices, box, Li_O_water_minima)
                if fraction == "No neighbors":
                    windows_counted-=1
                    start_frame+= skip_step
                    continue
                fractions.append(fraction)
                fractions_mean=((windows_counted-1)*fractions_mean+fractions[-1])/windows_counted
                start_frame+= skip_step

            fraction_overall_list.append(fractions_mean)

        fraction_overall = np.mean(fraction_overall_list, axis = 0).T
        t = np.array(list(range(len(fraction_overall))))*5/1000
        index = find_nearest(fraction_overall, 0.5)[0]
        print("Approx time to lose one neighbor is {}".format(t[index]))
         
        print("The size of overall is ", np.shape(fraction_overall))
        np.savetxt(
            f"fraction.txt",
            fraction_overall,
        )



        fig = plt.figure()
        ax = fig.gca()
        ax.scatter(t, fraction_overall, label="Fraction", s=0.5)
        #ax.plot(t, single_exp(t, tau),'orange', label="y=exp(-t/{:.1f})".format( tau ),  )
        if fraction_overall[-1]<0.5:
            plt.axvline(x=t[index], ymin=0, ymax=1,color = 'magenta', label='Tau = {} ps'.format(t[index]))
        plt.legend()
        plt.grid(alpha=0.2)
        plt.xlabel("Time (ps)")
        plt.ylabel("Fraction")
        #plt.xlim([0,1000])
        plt.ylim([0,1])
        plt.tight_layout()
        plt.savefig("fraction.pdf")
        plt.close()
        os.chdir("..")




if __name__ == "__main__":
    main()

