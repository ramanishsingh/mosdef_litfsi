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
    data_path = "rdf_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)
    
    top="../../10m_small.pdb"
    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T","Functional")):
        #if temp<300 or functional =="BLYP":
        #    continue
        print(temp,functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list=[]
        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            full_traj = md.load(job.fn("litfsi-pos-1.xyz"), top=top)
            #print(full_traj.n_frames)
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

        output_string=""
        for frame in range(comb_traj.n_frames):
            output_string+="{}\n\n".format(num_residues)

            for residue in comb_traj.topology.residues:
        
                atom_names=[]
                atom_indices=[]
                atom_masses=[]

                for atom in residue.atoms:
                    atom_names.append(atom.element.symbol)
                    atom_masses.append(atom.element.mass)
                    atom_indices.append(atom.index)

                points=comb_traj.xyz[frame][atom_indices]
                com=box.center_of_mass(points,masses=atom_masses)
                output_string+="{} ".format(residue.name)+str(10*com[0])+" "+str(10*com[1])+" "+str(10*com[2])+"\n"
            if frame==0:
                text_file = open("frame_for_mbuild.xyz", "w")
                n = text_file.write(output_string)
                text_file.close()

        text_file = open("com_coors.xyz", "w")
        n = text_file.write(output_string)
        text_file.close()

        mb_frame=mb.load('frame_for_mbuild.xyz')
        mb_frame.save('com_top.pdb')

        com_traj=md.load("com_coors.xyz", top='com_top.pdb')
        com_traj= md.Trajectory(
                com_traj.xyz,
                com_traj.top,
                unitcell_lengths = np.tile([length[0]/10, length[1]/10, length[2]/10], (comb_traj.n_frames,1)),
                unitcell_angles = np.tile([90.,90.,90.], (comb_traj.n_frames,1)),
            )

        # com_traj  is the COM traj and comb_traj is the complete traj, now find the rdfs
        r_min=0.01
        r_max=min(length)/20-0.01
        n_bins=150
        print("r_max is {} nm".format(r_max))
       
        # third rdf is TFSI-water (COM-COM)
        pair_indices_1=com_traj.top.select("name TF2")
        #print(pair_indices_1)
        pair_indices_2=com_traj.top.select("name wat")
        #print(pair_indices_2)

        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(com_traj.unitcell_vectors), com_traj.xyz[:, pair_indices_1, :],com_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        plt.grid(alpha=0.2)
        plt.savefig("TFSI-water.png")
        np.savetxt(
            f"TFSI-water.txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.savefig("cdf_TFSI-water.png")
        plt.close()
        # first rdf is C(TFSI)-C(TFSI)
        pair_indices_1=comb_traj.top.select("name C1 or name C2")
        #print(pair_indices_1)
        pair_indices_2=pair_indices_1
        
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        plt.grid(alpha=0.2)
        plt.savefig("C(TFSI)-C(TFSI).png")
        np.savetxt(
            f"C(TFSI)-C(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.savefig("cdf_C(TFSI)-C(TFSI).png")
        plt.close()

        # second rdf is TFSI-TFSI (COM-COM)
        pair_indices_1=com_traj.top.select("name TF2")
        #print(pair_indices_1)
        pair_indices_2=pair_indices_1

        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(com_traj.unitcell_vectors), com_traj.xyz[:, pair_indices_1, :],com_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        plt.grid(alpha=0.2)
        plt.savefig("TFSI-TFSI.png")
        
        np.savetxt(
            f"TFSI-TFSI.txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.savefig("cdf_TFSI-TFSI.png")
        plt.close()


        #fourth rdf is O(water)-O(water)
        pair_indices_1=comb_traj.top.select("resname wat and element O")
        #print(pair_indices_1)
        pair_indices_2=pair_indices_1

        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        plt.grid(alpha=0.2)
        plt.savefig("O(water)-O(water).png")
        np.savetxt(
            f"O(water)-O(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.savefig("cdf_O(water)-O(water).png")
        plt.close()


        #fifth rdf is Li-O(TFSI)
        pair_indices_1=comb_traj.top.select("element Li")
        #print(pair_indices_1)
        pair_indices_2=comb_traj.top.select("resname TF2 and element O")

        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        plt.grid(alpha=0.2)
        plt.savefig("Li-O(TFSI).png")
        np.savetxt(
            f"Li-O(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        Li_O_TFSI_minima=minima_in_range(r*10, gr, 2,4)[0]
        print("The first shell minima of Li-O(TFSI) lies at {} \AA".format(Li_O_TFSI_minima))
        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.savefig("cdf_Li-O(TFSI).png")
        plt.close()


        #sixth rdf is Li-O(water)
        pair_indices_1=comb_traj.top.select("element Li")
        #print(pair_indices_1)
        pair_indices_2=comb_traj.top.select("resname wat and element O")

        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        plt.grid(alpha=0.2)
        plt.savefig("Li-O(water).png")
        np.savetxt(
            f"Li-O(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        Li_O_water_minima=minima_in_range(r*10, gr, 2,4)[0]
        print("The first shell minima of Li-O(water) lies at {} \AA".format(Li_O_water_minima))
        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.savefig("cdf_Li-O(water).png")
        plt.close()

       #cluster analysis
        wat_O_indices=comb_traj.top.select("element O and resname wat")
        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            #print('The frame is {}/{}'.format(frame, comb_traj.n_frames))
            num_neighbors=[]
            points=comb_traj.xyz[frame][wat_O_indices]
            aq = freud.locality.AABBQuery(box, points)
            for point in range(points.shape[0]):
                #print('the point is {}'.format(point))
                neighbor_list=[]
                query_result = aq.query(points[point], dict(r_max=0.33))
                nlist = query_result.toNeighborList()
                query_result=0
                for (i,j) in nlist:
    #                print('nlist is {},{}'.format(i,j))
                    #neighbor_list.append(i)
                    neighbor_list.append(j)
                #print(neighbor_list)
                neighbor_list = [x for x in neighbor_list if x != point]
                #neighbor_list = [x for x in a if x != ]
   #             print('the neighbor list is {}'.format(neighbor_list))
                #print(neighbor_list)
                num_neighbors.append(len(neighbor_list))
                super_num_neighbors.append(len(neighbor_list))
        unique_elements=unique(super_num_neighbors)
        print(unique_elements)
        occur_freq=[]
        for element in unique_elements:
            occur_freq.append(super_num_neighbors.count(element))
        i=0
        for freq in occur_freq:
           occur_freq[i]=occur_freq[i]/len(super_num_neighbors)
           i+=1
        plt.figure()
        plt.scatter(unique_elements, occur_freq)
        plt.grid(alpha=0.2)
        unique_elements=np.array(unique_elements)
        occur_freq=np.array(occur_freq)
        occur_freq=np.reshape(occur_freq, (occur_freq.shape[0],1))
        unique_elements=np.reshape(unique_elements, (unique_elements.shape[0],1))
        c=np.concatenate((unique_elements,occur_freq),axis=1)
        ind=np.argsort(c[:,0])
        c=c[ind]
        plt.grid(alpha=0.2)
        np.savetxt('prob_O(water)_O(water)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_O(water)_O(water)_cluster.jpg')
        plt.close()

        wat_O_indices=comb_traj.top.select("element O and resname wat")
        Li_indices=comb_traj.top.select("element Li")

        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            #print('The frame is {}/{}'.format(frame, comb_traj.n_frames))
            num_neighbors=[]
            points_Li=comb_traj.xyz[frame][Li_indices]
            points_O=comb_traj.xyz[frame][wat_O_indices]
            aq = freud.locality.AABBQuery(box, points_O)
            for point in range(points_Li.shape[0]):
                #print('the point is {}'.format(point))
                neighbor_list=[]
                query_result = aq.query(points_Li[point], dict(r_max=Li_O_water_minima/10))
                nlist = query_result.toNeighborList()
                query_result=0
                for (i,j) in nlist:
    #                print('nlist is {},{}'.format(i,j))
                    #neighbor_list.append(i)
                    neighbor_list.append(j)
                #print(neighbor_list)
                neighbor_list = [x for x in neighbor_list if x != point]
                #neighbor_list = [x for x in a if x != ]
   #             print('the neighbor list is {}'.format(neighbor_list))
                #print(neighbor_list)
                num_neighbors.append(len(neighbor_list))
                super_num_neighbors.append(len(neighbor_list))
        unique_elements=unique(super_num_neighbors)
        print(unique_elements)
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
        np.savetxt('prob_Li_O(water)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_Li_O(water)_cluster.jpg')
        plt.close()

        TF2_O_indices=comb_traj.top.select("element O and resname TF2")
        Li_indices=comb_traj.top.select("element Li")

        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            #print('The frame is {}/{}'.format(frame, comb_traj.n_frames))
            num_neighbors=[]
            points_Li=comb_traj.xyz[frame][Li_indices]
            points_O=comb_traj.xyz[frame][TF2_O_indices]
            aq = freud.locality.AABBQuery(box, points_O)
            for point in range(points_Li.shape[0]):
                #print('the point is {}'.format(point))
                neighbor_list=[]
                query_result = aq.query(points_Li[point], dict(r_max=Li_O_TFSI_minima/10))
                nlist = query_result.toNeighborList()
                query_result=0
                for (i,j) in nlist:
    #                print('nlist is {},{}'.format(i,j))
                    #neighbor_list.append(i)
                    neighbor_list.append(j)
                #print(neighbor_list)
                neighbor_list = [x for x in neighbor_list if x != point]
                #neighbor_list = [x for x in a if x != ]
   #             print('the neighbor list is {}'.format(neighbor_list))
                #print(neighbor_list)
                num_neighbors.append(len(neighbor_list))
                super_num_neighbors.append(len(neighbor_list))
        unique_elements=unique(super_num_neighbors)
        print(unique_elements)
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
        np.savetxt('prob_Li_O(TF2)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_Li_O(TF2)_cluster.jpg')
        plt.close()


        os.chdir("..")    
if __name__ == "__main__":
    main()
