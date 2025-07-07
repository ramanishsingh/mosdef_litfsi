import os
import freud
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import matplotlib
import shutil
import signac

from hbond_elliptical import calualateHBMap

def main():
    # All output data will be stored in rdf_data folder
    save_string_file = "" # This string will be saved to a file
    data_path = "hbond_elliptical_data_tfsi"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    # top file to be used for the project
    top="../../10m_small.pdb"

    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T","Functional")):
        #Only computing RDFs for PBE functional
        if functional =="BLYP":
            continue
        print(temp,functional)
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list = []
        seed_list = []
        for job in group:
            seed=job.sp.Seed
            seed_list.append(seed)
            os.mkdir(str(seed))
            length=job.sp.L
            full_traj = md.load(job.fn("litfsi-pos-1.xyz"), top=top)
            full_traj = md.Trajectory(
                    full_traj.xyz,
                    full_traj.top,
                    unitcell_lengths = np.tile([length[0]/10, length[1]/10, length[2]/10], (full_traj.n_frames,1)),
                    unitcell_angles = np.tile([90.,90.,90.], (full_traj.n_frames,1)),
                )
            if full_traj.n_frames>5000:

                traj_list.append(full_traj[5000:]) # discarding first 25 ps

        box = freud.box.Box(Lx=length[0]/10,Ly= length[1]/10, Lz=length[2]/10)
        top = full_traj.topology


        nbins_r = 400
        nbins_a = 400
        r_cutoff = 0.75

        skip_every_x_frames = 10
        sel_oxygen_head = 'resname TF2 and element O' ; sel_hydrogen = 'name H1 or name H2' ; sel_oxygen_tail = 'name O5'; list_names_hydrogen = ["H1", "H2"] ; list_names_oxygen_head = ["O1", "O2", "O3", "O4"] ; list_names_oxygen_tail = ["O5"]

        
        num_hbond_list = []
        for index, traj in enumerate(traj_list):
            seed = seed_list[index]
            os.chdir("{}".format(seed))


            rdf_output, inter_output, map_output,hbond,hbond_time = calualateHBMap(traj, r_cutoff, nbins_r, nbins_a, skip_every_x_frames, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail, bonded_pdb_provided = False)
            inter_output[0] = inter_output[0]*180/np.pi
            print("The max element is ", np.amax(map_output))

            plt.figure()
            cmap = plt.get_cmap('jet')
            plt.figure(figsize=(5, 3))
            plt.style.use('default')
            levels = np.linspace(1,20,11)
            #levels = np.linspace(0,10,11) 
            cs = plt.contourf(rdf_output[0], inter_output[0], map_output,cmap=cmap, levels= levels)
            plt.xlabel('$r$ (nm)')
            plt.ylabel('$\u03B8$ (degrees)')
    
            
        #z = np.linspace(-1,-0.5,100)
        #x = np.linspace(0,0.4,100)
        #x,z = np.meshgrid(x,z)
        
        #ellipse = matplotlib.patches.Ellipse(
         #       (0.3, -1), 0.05*2, 0.3572*2, ec="magenta", facecolor="none", linewidth=3
         #   )
        #plt.gca().add_patch(ellipse)
        
        #f = ((x-0.30)/0.05)**2 + (1/0.3572**2)*(z+1)**2-1
        #plt.contour(x,z,f,[0, 0.01])
            
            plt.xlim([0.2, 0.4])
            plt.ylim([130, 180])
            plt.colorbar()
            plt.tight_layout()
            plt.savefig("hbond.png")                

            hbond_time_array  = np.array(hbond_time)
            np.savetxt("hbond_time.txt", hbond_time_array)
            save_string_file+="The avg number of hbonds per frame is {} \n".format(np.mean(hbond_time_array))
            np.savetxt('map_output.csv',map_output,delimiter=",")
            np.savetxt('r.csv', rdf_output[0], delimiter=",")
            np.savetxt('theta.csv', inter_output[0], delimiter=",")
            num_hbond_list.append(np.mean(hbond_time_array))
            os.chdir("..")

        os.chdir("..")
        save_string_file+= "overall hbond mean is {} and sem is {}".format(np.mean(num_hbond_list), np.std(num_hbond_list)/(len(num_hbond_list)**0.5))
        text_file = open("hbond_analysis_info_{}_{}.txt".format(temp, functional), "w")
        n = text_file.write(save_string_file)
        text_file.close()
        print(num_hbond_list)

if __name__ == "__main__":
    main()

