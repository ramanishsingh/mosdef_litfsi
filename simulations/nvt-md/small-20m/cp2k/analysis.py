import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import signac
import shutil
from scipy import stats
import freud
import mbuild as mb



def main():
    data_path = "rdf_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)
    
    top="../../20m_small.pdb"
    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T","Functional")):
        print(temp,functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list=[]
        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            full_traj = md.load(job.fn("litfsi-pos-1.xyz"), top=top)
            print(full_traj.n_frames)
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
        print(comb_traj)
        print("The combined trajectory has {} frames = {} ps ".format(comb_traj.n_frames,comb_traj.n_frames*5/1000))
        print('All residues: %s' % [residue for residue in comb_traj.topology.residues])
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
        # com_traj  is the COM traj and comb_traj is the complete traj, now find the rdfs
        r_min=0.01
        r_max=min(length)/20-0.01
        n_bins=300
        print("r_max is {} nm".format(r_max))
       
        # first rdf is C(TFSI)-C(TFSI)
        pair_indices_1=comb_traj.top.select("name C1 or name C2")
        print(pair_indices_1)
        pair_indices_2=pair_indices_1
        
        freud_rdf = freud.density.RDF(bins=n_bins, r_min=r_min, r_max=r_max)
        for a,b,c in zip(np.asarray(comb_traj.unitcell_vectors), comb_traj.xyz[:, pair_indices_1, :],comb_traj.xyz[:, pair_indices_2, :]):
            freud_rdf.compute((a,b),query_points=c, reset=False)
        r, gr, cdf = freud_rdf.bin_centers, freud_rdf.rdf, freud_rdf.n_r
        plt.figure()
        plt.plot(r*10,gr)
        plt.savefig("C(TFSI)-C(TFSI).png")
        plt.close()
             
        os.chdir("..")    
if __name__ == "__main__":
    main()
