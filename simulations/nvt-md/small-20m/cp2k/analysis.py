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

def snap_molecule_indices(system, r_max):
    """Find molecule index for each particle.

    Given a snapshot from a trajectory, compute clusters of bonded molecules
    and return an array of the molecule index of each particle.

    Parameters
    ----------
    snap : (box, points)

    Returns
    -------
    numpy array (N_particles,)

    """
    system = freud.AABBQuery.from_system(system)
    cl = freud.cluster.Cluster()
    cl.compute(system, neighbors={'r_max': r_max})
    

    return cl.cluster_idx



def main():
    # All output data will be stored in rdf_data folder
    save_string_file = "" # This string will be saved to a file 
    data_path = "rdf_data"
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
        comb_traj.save("comb_traj.xyz")
        #comb_traj.save("comb_traj.pdb")
        save_string_file += "\n"+"The combined trajectory has {} frames = {} ps ".format(comb_traj.n_frames,comb_traj.n_frames*5/1000)+"\n"
        box = freud.box.Box(Lx=length[0]/10,Ly= length[1]/10, Lz=length[2]/10)
        num_residues=comb_traj.n_residues
        
        output_string="" #This output_string contains the xyz file for the com_traj

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
       
        # first rdf is TFSI-water (COM-COM)
        pair_indices_1=com_traj.top.select("name TF2")
        pair_indices_2=com_traj.top.select("name wat")

        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(com_traj.unitcell_vectors), com_traj.xyz[:, pair_indices_1, :],com_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
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
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_TFSI-water.png")
        plt.close()


        # second rdf is C(TFSI)-C(TFSI)
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
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
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
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_C(TFSI)-C(TFSI).png")
        plt.close()

        # third rdf is TFSI-TFSI (COM-COM)
        pair_indices_1=com_traj.top.select("name TF2")
        pair_indices_2=pair_indices_1

        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(com_traj.unitcell_vectors), com_traj.xyz[:, pair_indices_1, :],com_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
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
        pair_indices_2=pair_indices_1

        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        O_O_minima=minima_in_range(r*10, gr, 3,4.5)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=O_O_minima, color='orange', label='Minima at {}'.format(round(O_O_minima,2)), ls='--')
        plt.figure()
        plt.plot(r*10,gr)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("O(water)-O(water).png")
        np.savetxt(
            f"O(water)-O(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of O-O lies at {} \AA".format(O_O_minima))
        save_string_file += "\n"+ "The first shell minima of O-O lies at {} \AA".format(O_O_minima) +"\n"


        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_O(water)-O(water).png")
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
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_Li-O(TFSI).png")
        plt.close()


        #sixth rdf is Li-O(water)
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
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_Li-N(TFSI).png")
        plt.close()


        #eighth rdf is Li-F(TFSI)
        pair_indices_1=comb_traj.top.select("element Li")
        pair_indices_2=comb_traj.top.select("resname TF2 and element F")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        Li_F_TFSI_minima=minima_in_range(r*10, gr, 2,4)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=Li_F_TFSI_minima, color='orange', label='Minima at {}'.format(round(Li_F_TFSI_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("Li-F(TFSI).png")
        np.savetxt(
            f"Li-F(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of Li-F(TFSI) lies at {} \AA".format(Li_F_TFSI_minima))
        save_string_file += "\n"+ "The first shell minima of Li-F(TFSI) lies at {} \AA".format(Li_F_TFSI_minima) +"\n"

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_Li-F(TFSI).png")
        plt.close()


        #ninth rdf is Li-Li
        pair_indices_1=comb_traj.top.select("element Li")
        pair_indices_2=comb_traj.top.select("element Li")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        Li_Li_minima=minima_in_range(r*10, gr, 3,5)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=Li_Li_minima, color='orange', label='Minima at {}'.format(round(Li_Li_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("Li-Li.png")
        np.savetxt(
            f"Li-Li.txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of Li-Li lies at {} \AA".format(Li_Li_minima))
        save_string_file += "\n"+ "The first shell minima of Li-Li lies at {} \AA".format(Li_Li_minima) +"\n"

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_Li-Li.png")
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
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_Li-S(TFSI).png")
        plt.close()

        #eleventh rdf is S(TFSI)-S(TFSI)
        pair_indices_1=comb_traj.top.select("resname TF2 and element S")
        pair_indices_2=comb_traj.top.select("resname TF2 and element S")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        S_S_TFSI_minima=minima_in_range(r*10, gr, 3.62,5)[0]
        #make a vertical line at the minima in the plot
        #plt.axvline(x=S_S_TFSI_minima, color='orange', label='Minima at {}'.format(round(S_S_TFSI_minima,2)), ls='--')
        #plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("S-S(TFSI).png")
        np.savetxt(
            f"S-S(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        #print("The first shell minima of S-S(TFSI) lies at {} \AA".format(S_S_TFSI_minima))
        #save_string_file += "\n"+ "The first shell minima of S-S(TFSI) lies at {} \AA".format(Li_S_TFSI_minima) +"\n"

        ##########

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_S(TFSI)-S(TFSI).png")
        plt.close()

        #twelfth RDF is O(TFSI)-H(water)

        pair_indices_1=comb_traj.top.select("resname TF2 and element O")
        pair_indices_2=comb_traj.top.select("resname wat and element H")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        O_TFSI_H_water_minima=minima_in_range(r*10, gr, 2,3)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=O_TFSI_H_water_minima, color='orange', label='Minima at {}'.format(round(O_TFSI_H_water_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("O(TFSI)-H(water).png")
        np.savetxt(
            f"O(TFSI)-H(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of O(TFSI)-H(water) lies at {} \AA".format(O_TFSI_H_water_minima))
        save_string_file += "\n"+ "The first shell minima of O(TFSI)-H(water) lies at {} \AA".format(O_TFSI_H_water_minima) +"\n"

        ##########

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_O(TFSI)-H(water).png")
        plt.close()


        #thirteenth RDF is O(TFSI)-O(water)

        pair_indices_1=comb_traj.top.select("resname TF2 and element O")
        pair_indices_2=comb_traj.top.select("resname wat and element O")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        O_TFSI_O_water_minima=minima_in_range(r*10, gr, 3.1,4.9)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=O_TFSI_O_water_minima, color='orange', label='Minima at {}'.format(round(O_TFSI_O_water_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("O(TFSI)-O(water).png")
        np.savetxt(
            f"O(TFSI)-O(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of O(TFSI)-O(water) lies at {} \AA".format(O_TFSI_O_water_minima))
        save_string_file += "\n"+ "The first shell minima of O(TFSI)-O(water) lies at {} \AA".format(O_TFSI_O_water_minima) +"\n"
        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_O(TFSI)-O(water).png")
        plt.close()



        #fourteenth RDF is F(TFSI)-H(water)

        pair_indices_1=comb_traj.top.select("resname TF2 and element F")
        pair_indices_2=comb_traj.top.select("resname wat and element H")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        F_TFSI_H_water_minima=minima_in_range(r*10, gr, 2,3)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=F_TFSI_H_water_minima, color='orange', label='Minima at {}'.format(round(F_TFSI_H_water_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("F(TFSI)-H(water).png")
        np.savetxt(
            f"F(TFSI)-H(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of F(TFSI)-H(water) lies at {} \AA".format(F_TFSI_H_water_minima))
        save_string_file += "\n"+ "The first shell minima of F(TFSI)-H(water) lies at {} \AA".format(F_TFSI_H_water_minima) +"\n"

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_F(TFSI)-H(water).png")
        plt.close()



        #fifteenth RDF is N(TFSI)-H(water)

        pair_indices_1=comb_traj.top.select("resname TF2 and element N")
        pair_indices_2=comb_traj.top.select("resname wat and element H")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        N_TFSI_H_water_minima=minima_in_range(r*10, gr, 2,3)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=N_TFSI_H_water_minima, color='orange', label='Minima at {}'.format(round(N_TFSI_H_water_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("N(TFSI)-H(water).png")
        np.savetxt(
            f"N(TFSI)-H(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of N(TFSI)-H(water) lies at {} \AA".format(N_TFSI_H_water_minima))
        save_string_file += "\n"+ "The first shell minima of N(TFSI)-H(water) lies at {} \AA".format(N_TFSI_H_water_minima) +"\n"
        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_N(TFSI)-H(water).png")
        plt.close()



        #sixteenth RDF is F(TFSI)-O(water)

        pair_indices_1=comb_traj.top.select("resname TF2 and element F")
        pair_indices_2=comb_traj.top.select("resname wat and element O")
        gr=0
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        F_TFSI_O_water_minima=minima_in_range(r*10, gr, 4,5)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=F_TFSI_O_water_minima, color='orange', label='Minima at {}'.format(round(F_TFSI_H_water_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("F(TFSI)_O(water).png")
        np.savetxt(
            f"F(TFSI)-O(water).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of F(TFSI)-O(water) lies at {} \AA".format(F_TFSI_O_water_minima))
        save_string_file += "\n"+ "The first shell minima of F(TFSI)-O(water) lies at {} \AA".format(F_TFSI_O_water_minima) +"\n"
        ##########

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_F(TFSI)-O(water).png")
        plt.close()



       ##seventeenth RDF is intermolecular S(TFSI) -S(TFSI)

        S_indices = comb_traj.top.select("resname TF2 and element S")
        S_traj = comb_traj.xyz[:, S_indices, :]

        system = (box, S_traj[0])

        molecules = snap_molecule_indices(system, r_max = 0.4)
        print("molecules found for S(TFSI) are ", molecules)
        freud_rdf = freud.density.RDF(bins=n_bins, r_max=r_max, r_min=r_min)
        for frame in range(len(S_traj)):
            system = (box, S_traj[frame])
            aq = freud.locality.AABBQuery.from_system(system)
            nlist = aq.query(
                S_traj[frame], {"r_max": r_max, "exclude_ii": True}
            ).toNeighborList()


            pre_filter = len(nlist)
            indices_A = molecules[nlist.point_indices]
            indices_B = molecules[nlist.query_point_indices]
            nlist.filter(indices_A != indices_B)
            post_filter = len(nlist)

            freud_rdf.compute(aq, neighbors=nlist, reset=False)
            normalization = post_filter / pre_filter 
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf*normalization, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        S_TFSI_S_TFSI_minima=minima_in_range(r*10, gr, 4,6)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=S_TFSI_S_TFSI_minima, color='orange', label='Minima at {}'.format(round(S_TFSI_S_TFSI_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("intermolecular-S(TFSI)-S(TFSI).png")
        np.savetxt(
            f"intermolecular-S(TFSI)-S(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of intermolecular-S(TFSI)-S(TFSI) lies at {} \AA".format(S_TFSI_S_TFSI_minima))
        save_string_file += "\n"+ "The first shell minima of S(TFSI)-S(TFSI) lies at {} \AA".format(S_TFSI_S_TFSI_minima) +"\n"
        ##########

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_intermolecular-S(TFSI)-S(TFSI).png")
        plt.close()

    

       ##eighteenth RDF is intermolecular C(TFSI) -C(TFSI)

        C_indices = comb_traj.top.select("resname TF2 and element C")
        C_traj = comb_traj.xyz[:, C_indices, :]

        system = (box, C_traj[0])

        #molecules = snap_molecule_indices(system, r_max = 0.4)
        print("molecules found for C(TFSI) are ", molecules)
        freud_rdf = freud.density.RDF(bins=n_bins, r_max=r_max, r_min=r_min)
        for frame in range(len(S_traj)):
            system = (box, C_traj[frame])
            aq = freud.locality.AABBQuery.from_system(system)
            nlist = aq.query(
                C_traj[frame], {"r_max": r_max, "exclude_ii": True}
            ).toNeighborList()


            pre_filter = len(nlist)
            indices_A = molecules[nlist.point_indices]
            indices_B = molecules[nlist.query_point_indices]
            nlist.filter(indices_A != indices_B)
            post_filter = len(nlist)

            freud_rdf.compute(aq, neighbors=nlist, reset=False)
            normalization = post_filter / pre_filter
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf*normalization, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        C_TFSI_C_TFSI_minima=minima_in_range(r*10, gr, 6,r_max*10)[0]
        #make a vertical line at the minima in the plot
        plt.axvline(x=C_TFSI_C_TFSI_minima, color='orange', label='Minima at {}'.format(round(C_TFSI_C_TFSI_minima,2)), ls='--')
        plt.legend()
        plt.grid(alpha=0.2)
        plt.legend()
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$g(r)$")
        plt.savefig("intermolecular-C(TFSI)-C(TFSI).png")
        np.savetxt(
            f"intermolecular-C(TFSI)-C(TFSI).txt",
            np.transpose(np.vstack([r*10, gr,r*10,cdf])),
            header="distance (\AA)\tRDF\tdistance (\AA)\tcdf",
        )
        print("The first shell minima of intermolecular-C(TFSI)-C(TFSI) lies at {} \AA".format(C_TFSI_C_TFSI_minima))
        save_string_file += "\n"+ "The first shell minima of intermolecular C(TFSI)-C(TFSI) lies at {} \AA".format(C_TFSI_C_TFSI_minima) +"\n"
        ##########

        plt.close()
        plt.figure()
        plt.plot(r*10,cdf)
        plt.grid(alpha=0.2)
        plt.xlabel(r"$r$/$\mathrm{\AA}$")
        plt.ylabel("$n(r)$")
        plt.savefig("cdf_intermolecular-C(TFSI)-C(TFSI).png")
        plt.close()

   
       #cluster analysis

        wat_O_indices=comb_traj.top.select("element O and resname wat")
        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            num_neighbors=[]
            points=comb_traj.xyz[frame][wat_O_indices]
            aq = freud.locality.AABBQuery(box, points)
            for point in range(points.shape[0]):
                #print('the point is {}'.format(point))
                neighbor_list=[]
                query_result = aq.query(points[point], dict(r_max=0.33, r_min = 0.01))
                nlist = query_result.toNeighborList()
                query_result=0
                for (i,j) in nlist:
                    neighbor_list.append(j)
                neighbor_list = [x for x in neighbor_list if x != point]
                num_neighbors.append(len(neighbor_list))
                super_num_neighbors.append(len(neighbor_list))
        unique_elements=unique(super_num_neighbors)
       # print(unique_elements)
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
        plt.xlabel("Number of Oxygens")
        plt.ylabel("Probability")
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
            num_neighbors=[]
            points_Li=comb_traj.xyz[frame][Li_indices]
            points_O=comb_traj.xyz[frame][wat_O_indices]
            aq = freud.locality.AABBQuery(box, points_O)
            for point in range(points_Li.shape[0]):
                neighbor_list=[]
                query_result = aq.query(points_Li[point], dict(r_max=Li_O_water_minima/10))
                nlist = query_result.toNeighborList()
                query_result=0
                #for (i,j) in nlist:
                #    neighbor_list.append(j)
                #neighbor_list = [x for x in neighbor_list if x != point]
                num_neighbors.append(len(nlist))
                super_num_neighbors.append(len(nlist))
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
        np.savetxt('prob_Li_O(water)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_Li_O(water)_cluster.jpg')
        plt.close()


        TF2_O_indices=comb_traj.top.select("element O and resname TF2")
        Li_indices=comb_traj.top.select("element Li")

        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            num_neighbors=[]
            points_Li=comb_traj.xyz[frame][Li_indices]
            points_O=comb_traj.xyz[frame][TF2_O_indices]
            aq = freud.locality.AABBQuery(box, points_O)
            for point in range(points_Li.shape[0]):
                neighbor_list=[]
                query_result = aq.query(points_Li[point], dict(r_max=Li_O_TFSI_minima/10))
                nlist = query_result.toNeighborList()
                query_result=0
                #for (i,j) in nlist:
                #    neighbor_list.append(j)
                #neighbor_list = [x for x in neighbor_list if x != point]
                num_neighbors.append(len(nlist))
                super_num_neighbors.append(len(nlist))
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
        np.savetxt('prob_Li_O(TF2)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_Li_O(TF2)_cluster.jpg')
        plt.close()

        # O (water) -O (TF2) cluster

        TF2_O_indices=comb_traj.top.select("element O and resname TF2")
        wat_O_indices=comb_traj.top.select("element O and resname wat")

        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            num_neighbors=[]
            points_wat_O=comb_traj.xyz[frame][wat_O_indices]
            points_O=comb_traj.xyz[frame][TF2_O_indices]
            aq = freud.locality.AABBQuery(box, points_O)
            for point in range(points_wat_O.shape[0]):
                neighbor_list=[]
                query_result = aq.query(points_wat_O[point], dict(r_max=O_TFSI_O_water_minima/10))
                nlist = query_result.toNeighborList()
                query_result=0
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
        plt.xlabel("Number of Oxygens (TFSI)")
        plt.ylabel("Probability")
        np.savetxt('prob_O(wat)_O(TF2)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_O(wat)_O(TF2)_cluster.jpg')
        plt.close()



        # O (TF2) -O (wat) cluster

        TF2_O_indices=comb_traj.top.select("element O and resname TF2")
        wat_O_indices=comb_traj.top.select("element O and resname wat")

        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            num_neighbors=[]
            points_wat_O=comb_traj.xyz[frame][wat_O_indices]
            points_O=comb_traj.xyz[frame][TF2_O_indices]
            aq = freud.locality.AABBQuery(box, points_wat_O)
            for point in range(points_O.shape[0]):
                neighbor_list=[]
                query_result = aq.query(points_O[point], dict(r_max=O_TFSI_O_water_minima/10))
                nlist = query_result.toNeighborList()
                query_result=0
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
        plt.xlabel("Number of Oxygens (water)")
        plt.ylabel("Probability")
        np.savetxt('prob_O(TF2)-O(wat)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_O(TF2)-O(wat)_cluster.jpg')
        plt.close()



        # cluster size cumulative probability
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
        #print("The sum of all freq for one frame is {}".format(np.sum(frequency)/comb_traj.n_frames))
        frac_cluster_size=np.array(list(counter.keys()))
        p = frac_cluster_size.argsort()
        frac_cluster_size=frac_cluster_size[p]/points_O.shape[0]
        frequency=frequency[p]
        #print(frequency)
        prob=[]
        for i in range(frequency.shape[0]):
         #   print("i is {}".format(i))
         #   print("local frq is {}".format(frequency[i]))
         #   print("sum of all freq is {}".format(np.sum(frequency)))
         #   print("ratio is {}".format(frequency[i]/np.sum(frequency)))
            prob.append(frequency[i]*frac_cluster_size[i]/comb_traj.n_frames)
        prob=np.array(prob)
        cum_prob=np.cumsum(prob)
        #print("prob is  {} ".format(prob))
        #print("cum prob is  {} ".format(cum_prob))
        #print("frac cluster size is {}".format(frac_cluster_size))
        plt.figure()
        plt.loglog(frac_cluster_size,cum_prob)
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


        S_indices = comb_traj.top.select("resname TF2 and element S")
        S_traj = comb_traj.xyz[:, S_indices, :]
        system = (box, S_traj[0])
        molecules = snap_molecule_indices(system, r_max = 0.4)

        TFSI_S_indices=comb_traj.top.select("element S and resname TF2")

        super_cluster_size=[]
        super_num_clusters=[]
        for frame in range(comb_traj.n_frames):
        #for frame in range(1):
            num_neighbors=[]
            points_S=comb_traj.xyz[frame][TFSI_S_indices]
            points_S=box.wrap(points_S)
            system = freud.AABBQuery(box, points_S)
            cl = freud.cluster.Cluster()
            cl.compute(system, neighbors={"r_max": 0.56})
            atomic_clusters = []
            for cid in range(cl.num_clusters):
                atomic_clusters.append(cl.cluster_keys[cid])
            TFSI_clusters = []
            for atomic_cluster in atomic_clusters:
                TFSI_clusters.append(unique(molecules[atomic_cluster]))
            for TFSI_cluster in TFSI_clusters:
                super_cluster_size.append(len(TFSI_cluster))



        counter=collections.Counter(super_cluster_size)
        frequency_dict=dict(counter)
        frequency= np.array(list(counter.values()))
        #print("The sum of all freq for one frame is {}".format(np.sum(frequency)/comb_traj.n_frames))
        frac_cluster_size=np.array(list(counter.keys()))
        p = frac_cluster_size.argsort()
        frac_cluster_size=frac_cluster_size[p]/(points_S.shape[0]/2)
        frequency=frequency[p]
        #print(frequency)
        prob=[]
        for i in range(frequency.shape[0]):
            prob.append(frequency[i]*frac_cluster_size[i]/comb_traj.n_frames)
        prob=np.array(prob)
        cum_prob=np.cumsum(prob)
        #print("prob is  {} ".format(prob))
        #print("cum prob is  {} ".format(cum_prob))
        #print("frac cluster size is {}".format(frac_cluster_size))
        plt.figure()
        plt.loglog(frac_cluster_size,cum_prob)
        plt.grid(alpha=0.2)
        plt.xlim([0.0001, 1])
        plt.ylim([1e-6, 1])
        plt.xlabel("N_{TFSI in cluster}/N_{total TFSI in box}")
        plt.ylabel("Cumulative probability")

        plt.savefig("cum_prob_S_S.png")
        np.savetxt(
            f"cum_prob_S_S.txt",
            np.transpose(np.vstack([frac_cluster_size, cum_prob])),
            header="frac_cluster_size\tcum_prob",
        )

        plt.close()

        #TFSI cluster based on C


        C_indices = comb_traj.top.select("resname TF2 and element C")
        C_traj = comb_traj.xyz[:, C_indices, :]
        system = (box, C_traj[0])
        molecules = snap_molecule_indices(system, r_max = 0.4)

        TFSI_C_indices=comb_traj.top.select("element C and resname TF2")

        super_cluster_size=[]
        super_num_clusters=[]
        for frame in range(comb_traj.n_frames):
        #for frame in range(1):
            num_neighbors=[]
            points_C=comb_traj.xyz[frame][TFSI_C_indices]
            points_C=box.wrap(points_C)
            system = freud.AABBQuery(box, points_C)
            cl = freud.cluster.Cluster()
            cl.compute(system, neighbors={"r_max": 0.56})
            atomic_clusters = []
            for cid in range(cl.num_clusters):
                atomic_clusters.append(cl.cluster_keys[cid])
            TFSI_clusters = []
            for atomic_cluster in atomic_clusters:
                TFSI_clusters.append(unique(molecules[atomic_cluster]))
            for TFSI_cluster in TFSI_clusters:
                super_cluster_size.append(len(TFSI_cluster))



        counter=collections.Counter(super_cluster_size)
        frequency_dict=dict(counter)
        frequency= np.array(list(counter.values()))
        #print("The sum of all freq for one frame is {}".format(np.sum(frequency)/comb_traj.n_frames))
        frac_cluster_size=np.array(list(counter.keys()))
        p = frac_cluster_size.argsort()
        frac_cluster_size=frac_cluster_size[p]/(points_C.shape[0]/2)
        frequency=frequency[p]
        #print(frequency)
        prob=[]
        for i in range(frequency.shape[0]):
            prob.append(frequency[i]*frac_cluster_size[i]/comb_traj.n_frames)
        prob=np.array(prob)
        cum_prob=np.cumsum(prob)
        #print("prob is  {} ".format(prob))
        #print("cum prob is  {} ".format(cum_prob))
        #print("frac cluster size is {}".format(frac_cluster_size))
        plt.figure()
        plt.loglog(frac_cluster_size,cum_prob)
        plt.grid(alpha=0.2)
        plt.xlim([0.0001, 1])
        plt.ylim([1e-6, 1])
        plt.xlabel("N_{TFSI in cluster}/N_{total TFSI in box}")
        plt.ylabel("Cumulative probability")

        plt.savefig("cum_prob_C_C.png")
        np.savetxt(
            f"cum_prob_C_C.txt",
            np.transpose(np.vstack([frac_cluster_size, cum_prob])),
            header="frac_cluster_size\tcum_prob",
        )

        plt.close()


        os.chdir("..")    

        text_file = open("analysis_info_{}_{}.txt".format(temp, functional), "w")
        n = text_file.write(save_string_file)
        text_file.close()

if __name__ == "__main__":
    main()
