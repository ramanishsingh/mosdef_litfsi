import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import signac
import shutil
from scipy import stats
import freud
import mbuild as mb
from collections import Counter



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
    os.chdir(data_path)
    
    top="../../10m_small.pdb"
    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T","Functional")):
        print(temp,functional)
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


        wat_O_indices=comb_traj.top.select("element O and resname wat")
        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            print('The frame is {}/{}'.format(frame, comb_traj.n_frames))
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
        unique_elements=np.array(unique_elements)
        occur_freq=np.array(occur_freq)
        occur_freq=np.reshape(occur_freq, (occur_freq.shape[0],1))
        unique_elements=np.reshape(unique_elements, (unique_elements.shape[0],1))
        c=np.concatenate((unique_elements,occur_freq),axis=1)
        ind=np.argsort(c[:,0])
        c=c[ind]
        plt.grid()
        np.savetxt('prob_O(water)_O(water)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_O(water)_O(water)_cluster.jpg')
        plt.close()

        wat_O_indices=comb_traj.top.select("element O and resname wat")
        Li_indices=comb_traj.top.select("element Li")

        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            print('The frame is {}/{}'.format(frame, comb_traj.n_frames))
            num_neighbors=[]
            points_Li=comb_traj.xyz[frame][Li_indices]
            points_O=comb_traj.xyz[frame][wat_O_indices]
            aq = freud.locality.AABBQuery(box, points_O)
            for point in range(points_Li.shape[0]):
                #print('the point is {}'.format(point))
                neighbor_list=[]
                query_result = aq.query(points_Li[point], dict(r_max=0.25))
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
        plt.grid()
        np.savetxt('prob_Li_O(water)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_Li_O(water)_cluster.jpg')
        plt.close()

        TF2_O_indices=comb_traj.top.select("element O and resname TF2")
        Li_indices=comb_traj.top.select("element Li")

        super_num_neighbors=[]
        for frame in range(comb_traj.n_frames):
            print('The frame is {}/{}'.format(frame, comb_traj.n_frames))
            num_neighbors=[]
            points_Li=comb_traj.xyz[frame][Li_indices]
            points_O=comb_traj.xyz[frame][TF2_O_indices]
            aq = freud.locality.AABBQuery(box, points_O)
            for point in range(points_Li.shape[0]):
                #print('the point is {}'.format(point))
                neighbor_list=[]
                query_result = aq.query(points_Li[point], dict(r_max=0.25))
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
        plt.grid()
        np.savetxt('prob_Li_O(TF2)_cluster.txt',c,header="cluster_size   prob")
        plt.savefig('prob_Li_O(TF2)_cluster.jpg')
        plt.close()


        os.chdir("..")    
if __name__ == "__main__":
    main()
