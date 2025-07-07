import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import signac
import shutil
from scipy import stats
import freud
import mbuild as mb
from scattering.utils.features import minima_in_range
import collections
import json


def main():
    # All output data will be stored in rdf_data folder
    save_string_file = "" # This string will be saved to a file 
    data_path = "neighbor_dist"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)


    # top file to be used for the project
    top="../../../10m_small.pdb"
    project = signac.get_project()
    for (temp), group in project.groupby(("T")):
        #Only computing VHFs for PBE functional
        functional = "MD"
        print(temp,functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list=[]


        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            trj_file = os.path.join(job.workspace(), "sample.xtc")
            full_traj = md.load(trj_file, top=top)
            print("full traj has n_frames", full_traj.n_frames)
            if full_traj.n_frames>1:

                traj_list.append(full_traj) # discarding first 10 ns
                print(" traj has n_frames", full_traj.n_frames)
        #get the shortest traj and reduce all trajs to that length
        print("All trajs loaded")
        box = freud.box.Box(Lx=length[0]/10,Ly= length[1]/10, Lz=length[2]/10)




        comb_traj=md.join(traj_list)
        # Add unit cell information
        comb_traj = md.Trajectory(
                comb_traj.xyz,
                comb_traj.top,
                unitcell_lengths = np.tile([length[0]/10, length[1]/10, length[2]/10], (comb_traj.n_frames,1)),
                unitcell_angles = np.tile([90.,90.,90.], (comb_traj.n_frames,1)),
            )
        print("The combined trajectory has {} frames = {} ps ".format(comb_traj.n_frames,comb_traj.n_frames*5/1000))
        #comb_traj.save("comb_traj.xyz")
        #comb_traj.save("comb_traj.pdb")
        save_string_file += "\n"+"The combined trajectory has {} frames = {} ps ".format(comb_traj.n_frames,comb_traj.n_frames*5/1000)+"\n"
        box = freud.box.Box(Lx=length[0]/10,Ly= length[1]/10, Lz=length[2]/10)
        num_residues=comb_traj.n_residues
        
        output_string="" #This output_string contains the xyz file for the com_traj

        #sixth rdf is Li-O(water)
        #Using this to find the minima


        # com_traj  is the COM traj and comb_traj is the complete traj, now find the rdfs
        r_min=0.01
        r_max=min(length)/20-0.01
        n_bins=150
        print("r_max is {} nm".format(r_max))

        pair_indices_1=comb_traj.top.select("element Li")
        pair_indices_2=comb_traj.top.select("resname wat and element O")
        gr = 0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        Li_O_water_minima=minima_in_range(r*10, gr, 2,4)[0]
        plt.axvline(x=Li_O_water_minima, color='orange', label='Minima at {}'.format(round(Li_O_water_minima,2)), ls='--')
        plt.legend()

        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("Li-O(water).png")
        np.savetxt(
            f"Li-O(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of Li-O(water) lies at {} \AA".format(Li_O_water_minima))
        save_string_file += "\n"+"The first shell minima of Li-O(water) lies at {} \AA".format(Li_O_water_minima) +"\n"
        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_Li-O(water).png")
        plt.close()



        #fifth rdf is Li-O(TFSI)
        pair_indices_1=comb_traj.top.select("element Li")
        pair_indices_2=comb_traj.top.select("resname TF2 and element O")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        Li_O_TFSI_minima=minima_in_range(r*10, gr, 2,4)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=Li_O_TFSI_minima, color='orange', label='Minima at {}'.format(round(Li_O_TFSI_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("Li-O(TFSI).png")
        np.savetxt(
            f"Li-O(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of Li-O(TFSI) lies at {} \AA".format(Li_O_TFSI_minima))
        save_string_file += "\n"+ "The first shell minima of Li-O(TFSI) lies at {} \AA".format(Li_O_TFSI_minima) +"\n"
        plt.close()

        #tenth rdf is Li-S(TFSI)
        pair_indices_1=comb_traj.top.select("element Li")
        pair_indices_2=comb_traj.top.select("resname TF2 and element S")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        Li_S_TFSI_minima=minima_in_range(r*10, gr, 3,5)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=Li_S_TFSI_minima, color='orange', label='Minima at {}'.format(round(Li_S_TFSI_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("Li-S(TFSI).png")
        np.savetxt(
            f"Li-S(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of Li-S(TFSI) lies at {} \AA".format(Li_S_TFSI_minima))
        save_string_file += "\n"+ "The first shell minima of Li-S(TFSI) lies at {} \AA".format(Li_S_TFSI_minima) +"\n"
        plt.close()



        counts = {}

        for i in range(8):
            for j in range(8):
                counts[(i,j)] =0

        probs = {}
        for i in range(8):
            for j in range(8):
                probs[(i,j)] =0


        wat_O_indices=comb_traj.top.select("element O and resname wat")
        Li_indices=comb_traj.top.select("element Li")
        TFSI_O_indices = comb_traj.top.select("resname TF2 and element O")
        TFSI_S_indices = comb_traj.top.select("resname TF2 and element S")
        super_num_neighbors=[]


        total_frames = len(comb_traj)
        nblocks = 4
        blocksize = int(total_frames/nblocks)
        block_num=0
        while block_num < nblocks:

            local_traj = comb_traj[blocksize*block_num:blocksize*(block_num+1)]
            skip = 1
            for frame in range(0, local_traj.n_frames, skip):
                num_neighbors=[]
                points_Li=local_traj.xyz[frame][Li_indices]
                points_O=local_traj.xyz[frame][wat_O_indices]
                points_O_TFSI = local_traj.xyz[frame][TFSI_O_indices]
                points_S_TFSI = local_traj.xyz[frame][TFSI_S_indices]
                aq = freud.locality.AABBQuery(box, points_O)
                for point in range(points_Li.shape[0]):
                    neighbor_list=[]
                    query_result = aq.query(points_Li[point], dict(r_max=Li_O_water_minima/10))
                    nlist = query_result.toNeighborList()
    
    
                    aq1 = freud.locality.AABBQuery(box, points_O_TFSI)
                    query1_result = aq1.query(points_Li[point], dict(r_max=Li_O_TFSI_minima/10))
                    nlist1 = query1_result.toNeighborList()
                    counts[(len(nlist), len(nlist1))]+=1
    
            total_sum = np.sum(list(counts.values()))
            for key in counts.keys():
    
                probs[key] = counts[key]/total_sum
    
            for key in counts.keys():
                if probs[key] < 0.01:
                    del probs[key]
    
    
            print(probs)
            # Serialize data into file:
            file = open('neighbor_frac_{}.txt'.format(block_num), 'w')
            file.write(str(probs))
            file.close()
            block_num+=1


        os.chdir("..")    

        text_file = open("analysis_info_{}_{}.txt".format(temp, functional), "w")
        n = text_file.write(save_string_file)
        text_file.close()

if __name__ == "__main__":
    main()
