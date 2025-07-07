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

def unique(list1):

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list



def give_two_closest(unique_TFSI_molecules, top, comb_traj, frame, original_Li_index):
    #print("Unique Tfsi molecules", unique_TFSI_molecules)
    overall_min_distances = []
    overall_min_indices = []
    for unique_TFSI_molecule in unique_TFSI_molecules:
        for residue in top.residues:
            if residue.index == unique_TFSI_molecule:
                atoms = list(residue.atoms)
                atom_indices = []
                for atom in atoms:
                    if (atom.element.symbol == "O" or atom.element.symbol == "S"):
                        atom_indices.append(atom.index)
                distances = []
                for index in atom_indices:
                    distances.append(md.compute_distances(comb_traj[frame], np.array([[index, original_Li_index]]))[0][0])
                #print(distances)
                min_distance = min(distances)
                min_index = distances.index(min_distance)
                overall_min_distances.append(min_distance)
                overall_min_indices.append(atom_indices[min_index])
    #print(overall_min_indices)
    order = np.argsort(overall_min_indices)[:2]

    indices = [overall_min_indices[order[0]], overall_min_indices[order[1]]]
    molecules = [unique_TFSI_molecules[order[0]], unique_TFSI_molecules[order[1]]]
    elements = [list(top.atoms)[indices[0]].element.symbol, list(top.atoms)[indices[1]].element.symbol]

    return indices, molecules, elements


def main():
    # All output data will be stored in rdf_data folder
    save_string_file = "" # This string will be saved to a file 
    data_path = "2dhist_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    # top file to be used for the project    
    top="../../../20m_small.pdb"
    project = signac.get_project()
    for (temp), group in project.groupby(("T")):
        functional = "FF-UND"
        print(temp,functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list=[]



        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            full_traj = md.load(job.fn("mol.lammpstrj"), top=top)
            print("full traj has n_frames", full_traj.n_frames)
            if full_traj.n_frames>10000:
                full_traj = md.Trajectory(
                        full_traj.xyz,
                        full_traj.top,
                        unitcell_lengths = np.tile([length[0]/10, length[1]/10, length[2]/10], (full_traj.n_frames,1)),
                        unitcell_angles = np.tile([90.,90.,90.], (full_traj.n_frames,1)),
                    )


                traj_list.append(full_traj[10001:]) # discarding first 10 ns
        #get the shortest traj and reduce all trajs to that length
        print("All trajs loaded")

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


        #seventh rdf is Li-N(TFSI)
        pair_indices_1=comb_traj.top.select("element Li")
        pair_indices_2=comb_traj.top.select("resname TF2 and element N")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        Li_N_TFSI_minima=minima_in_range(r*10, gr, 2,4)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=Li_N_TFSI_minima, color='orange', label='Minima at {}'.format(round(Li_N_TFSI_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("Li-N(TFSI).png")
        np.savetxt(
            f"Li-N(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of Li-N(TFSI) lies at {} \AA".format(Li_N_TFSI_minima))
        save_string_file += "\n"+ "The first shell minima of Li-N(TFSI) lies at {} \AA".format(Li_N_TFSI_minima) +"\n"
        plt.close()








        wat_O_indices=comb_traj.top.select("element O and resname wat")
        Li_indices=comb_traj.top.select("element Li")
        TFSI_O_indices = comb_traj.top.select("resname TF2 and element O")
        TFSI_S_indices = comb_traj.top.select("resname TF2 and element S")
        super_num_neighbors=[]

        num_neighbors_of_neighbors = []
        rrcos_triple = []
        #for frame in range(comb_traj.n_frames):
        skip = 1
        q_overall = []
        for frame in range(0, comb_traj.n_frames, skip):
            num_neighbors=[]
            points_Li=comb_traj.xyz[frame][Li_indices]
            points_O=comb_traj.xyz[frame][wat_O_indices]
            points_O_TFSI = comb_traj.xyz[frame][TFSI_O_indices]
            points_S_TFSI = comb_traj.xyz[frame][TFSI_S_indices]
            aq = freud.locality.AABBQuery(box, points_O)
            for point in range(points_Li.shape[0]):
                neighbor_list=[]
                query_result = aq.query(points_Li[point], dict(r_max=Li_O_water_minima/10))
                nlist = query_result.toNeighborList()
                if len(list(nlist)) ==2 :
                    two_neighbors = []
                    for (i,j) in nlist:
                        two_neighbors.append(j)
                    original_wat_indices = wat_O_indices[two_neighbors]
                    original_Li_index = Li_indices[point]

                    normal_angle = md.compute_angles(comb_traj[frame], np.array([[original_wat_indices[0], original_Li_index, original_wat_indices[1]]]))[0][0]
                    cosangle = np.cos(normal_angle)
                    r1 = md.compute_distances(comb_traj[frame], np.array([[original_wat_indices[0], original_Li_index]]))[0][0]
                    r2 = md.compute_distances(comb_traj[frame], np.array([[original_wat_indices[1], original_Li_index]]))[0][0]
                    rrcos_triple.append([r1, r2, cosangle])

                    #Finding if that Li with 2 water neighbors has any TFSI neighbors?

                    #oxygen of TFSI
                    aq1 = freud.locality.AABBQuery(box, points_O_TFSI)
                    query1_result = aq1.query(points_Li[point], dict(r_max=Li_O_TFSI_minima/10))
                    nlist1 = query1_result.toNeighborList()
                    print("Li with 2 water neighbors has {} Oxygens from TFSI".format(len(nlist1)))
                    O_TFSI_neighbors = []
                    for (i,j) in nlist1:
                        O_TFSI_neighbors.append(j)
                    original_O_TFSI_neighbor_indices = TFSI_O_indices[O_TFSI_neighbors]
                    TFSI_molecules = []
                    for index in original_O_TFSI_neighbor_indices:
                        TFSI_molecules.append(list(comb_traj.topology.atoms)[index].residue.index)

                    #Sulphur of TFSI
                    aq2 = freud.locality.AABBQuery(box, points_S_TFSI)
                    query2_result = aq2.query(points_Li[point], dict(r_max=Li_S_TFSI_minima/10))
                    nlist2 = query2_result.toNeighborList()
                    print("Li with 2 water neighbors has {} Sulphurs from TFSI".format(len(nlist2)))

                    S_TFSI_neighbors = []
                    for (i,j) in nlist2:
                        S_TFSI_neighbors.append(j)
                    original_S_TFSI_neighbor_indices = TFSI_S_indices[S_TFSI_neighbors]

                    for index in original_S_TFSI_neighbor_indices:
                        TFSI_molecules.append(list(comb_traj.topology.atoms)[index].residue.index)

                    #print(TFSI_molecules)
                    #print(unique(TFSI_molecules), len(unique(TFSI_molecules)))
                    num_neighbors_of_neighbors.append(len(unique(TFSI_molecules)))
                    if len(unique(TFSI_molecules))>1:
                        neighbor_indices, neighbor_molecules, neighbor_atoms = give_two_closest(unique(TFSI_molecules), comb_traj.top, comb_traj, frame, original_Li_index)

                    ##The elements of the the tetrahedron
                    #original_Li_index at the center
                    #original_wat_indices[0] and [1] are the original two neighbors
                    #neighbor_indices[0] and neighbor_indices[1] are the new two neighbors

                    cosangles = []
                    normal_angle = md.compute_angles(comb_traj[frame], np.array([[original_wat_indices[0], original_Li_index, original_wat_indices[1]]]))[0][0]
                    cosangles.append(np.cos(normal_angle))
                    normal_angle = md.compute_angles(comb_traj[frame], np.array([[original_wat_indices[0], original_Li_index, neighbor_indices[0]]]))[0][0]
                    cosangles.append(np.cos(normal_angle))
                    normal_angle = md.compute_angles(comb_traj[frame], np.array([[original_wat_indices[0], original_Li_index, neighbor_indices[1]]]))[0][0]
                    cosangles.append(np.cos(normal_angle))
                    normal_angle = md.compute_angles(comb_traj[frame], np.array([[original_wat_indices[1], original_Li_index, neighbor_indices[0]]]))[0][0]
                    cosangles.append(np.cos(normal_angle))
                    normal_angle = md.compute_angles(comb_traj[frame], np.array([[original_wat_indices[1], original_Li_index, neighbor_indices[1]]]))[0][0]
                    cosangles.append(np.cos(normal_angle))
                    normal_angle = md.compute_angles(comb_traj[frame], np.array([[neighbor_indices[0], original_Li_index, neighbor_indices[1]]]))[0][0]
                    cosangles.append(np.cos(normal_angle))

                    q = 1 - (3/8)*((cosangles[0]+1/3)**2+(cosangles[1]+1/3)**2+(cosangles[2]+1/3)**2+(cosangles[3]+1/3)**2+(cosangles[4]+1/3)**2+(cosangles[5]+1/3)**2)
                    q_overall.append(q)

                    print("=======================================================================")

                for (i,j) in nlist:
                    neighbor_list.append(j)
                neighbor_list = [x for x in neighbor_list if x != point]
                num_neighbors.append(len(neighbor_list))
                super_num_neighbors.append(len(neighbor_list))
        unique_elements=unique(super_num_neighbors)
        #print(unique_elements)
        occur_freq=[]
        for element in unique_elements:
            occur_freq.append(super_num_neighbors.count(element))
        i=0
        for freq in occur_freq:
            occur_freq[i]=occur_freq[i]/len(super_num_neighbors)
            i+=1
        plt.figure()
        plt.scatter(unique_elements, occur_freq)
        unique_elements=np.array(unique_elements)
        occur_freq=np.array(occur_freq)
        occur_freq=np.reshape(occur_freq, (occur_freq.shape[0],1))
        unique_elements=np.reshape(unique_elements, (unique_elements.shape[0],1))
        c=np.concatenate((unique_elements,occur_freq),axis=1)
        ind=np.argsort(c[:,0])
        c=c[ind]
        plt.grid(alpha=0.2)
        plt.xlabel("Number of Oxygens")
        plt.ylabel("Probability")
        rrcos_triple = np.asarray(rrcos_triple)


        nbins_r = 40
        nbins_a = 40
        r_cutoff = np.round(Li_O_water_minima/10,2)
        r_cutoff_min = 0.15
        radii, dr = np.linspace(r_cutoff_min, r_cutoff, nbins_r+1, retstep=True)
        r_centers = 0.5 * (radii[1:] + radii[:-1])

        cosangles, dcosa = np.linspace(-1, 1, nbins_a+1, retstep=True)
        a_centers = 0.5 * (cosangles[1:] + cosangles[:-1])
        first = rrcos_triple[:,::2]
        second = rrcos_triple[:, 1:]
        double = np.vstack((first, second))

        maps = np.zeros((nbins_a,nbins_r))

        for sample in double:
            dist = sample[0]
            cosangle = sample[1]
            index_d = np.digitize(dist, radii)-1
            index_i = np.digitize(cosangle,cosangles)-1
            maps[index_i, index_d]+=1

        for i in range(0, nbins_a):


            for j in range(0, nbins_r):
                rho_r_theta = (r_centers[j]**2)*dr

                maps[i,j] = maps[i,j] / (rho_r_theta)

        maps = maps/(sum(sum(maps))*dcosa*dr)


        plt.figure() ; cmap = plt.get_cmap('jet') ; plt.figure(figsize=(5, 3)) ;
        plt.style.use('default'); #levels = np.linspace(1,55,10) ; 
        #cs = plt.contourf(r_centers, a_centers, maps,levels = levels, cmap=cmap) ; 
        cs = plt.contourf(r_centers, a_centers, maps, cmap=cmap) 
        plt.xlabel('$r$ (nm)') ; plt.ylabel('cos ($\u03B8$)') ; 
        plt.xlim([r_cutoff_min, r_cutoff]) ; plt.ylim([-1, 1]) ;
        plt.colorbar()
        plt.tight_layout()
        plt.savefig("2dhist-rtheta.png", dpi = 600)

        np.savetxt('r-theta-samples.txt',double,header="r (nm)   cos(theta)")
        np.savetxt('2dhist-map.txt',maps)
        np.savetxt('r.txt',r_centers,header="r center (nm)")
        np.savetxt('costheta.txt',a_centers,header="cosangle center ")

        plt.figure()
        q_overall = np.array(q_overall)
        a = plt.hist(q_overall, density=True, range=[0,1], bins = 100)
        centers = 0.5*(a[1][1:]+a[1][:-1])
        hist_info = np.vstack((centers, a[0]))
        plt.figure()
        plt.plot(hist_info[0], hist_info[1])
        np.savetxt("q4.txt", hist_info.T)
        plt.xlabel("$q_4$")
        plt.ylabel("Normalized probability")
        plt.savefig("q_4_neighbors.png")
        plt.close()
        plt.figure()
        plt.hist(num_neighbors_of_neighbors, density=True,)
        labels, counts = np.unique(num_neighbors_of_neighbors, return_counts=True)
        counts = counts/sum(counts)
        hist_info = np.vstack((labels, counts)).T
        np.savetxt("non_water_neighbors.txt", hist_info)

        plt.xlabel("Number of non water neighbors")
        plt.ylabel("Normalized probability")
        plt.savefig("num_non_water_neighbors.png")
        plt.close()
        os.chdir("..")    

        text_file = open("analysis_info_{}_{}.txt".format(temp, functional), "w")
        n = text_file.write(save_string_file)
        text_file.close()

if __name__ == "__main__":
    main()
