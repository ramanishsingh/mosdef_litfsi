import main_files.analysis as analysis
import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv as csv


# Note: You may want to use VMD to convert the GOMC generated pdb file into a single pdb frame (frame 0 only)
# and an xyz file containing all the frames
# This will speed up analsyis and not consume >100GB memory. You will need to use only the step 0 frame from VMD
# as the pdb file and all frames of the pdb file converted into an xyz file.  The new files will need to be named
# SPCE_PORE_NVT_10_BOX_0_water_only.xyz and SPCE_PORE_NVT_10_BOX_0_at_0M_steps_water_only.pdb for the xyz and pdb file
# respectively.  This will need to be done for all 5 runs (i.e., 1r1, 1r2, 1r3, 1r4, 1r5), and changed as in the
# use_pdb_only_or_pdb_and_xyz varabile below.
# indicated in the notes below

use_pdb_only_or_pdb_and_xyz = 'pdb_only' # choose 'pdb_only' or 'pdb_and_xyz'

filepath_list = ['../set1/1r1', '../set1/1r2', '../set1/1r3', '../set1/1r4', '../set1/1r5']  #without slash after

pore_center = (5 + 5) / 10


for i in range(len(filepath_list)):
    filepath = filepath_list[i]

    print('Running ' + str(filepath_list[i]) + ' data:')

    if use_pdb_only_or_pdb_and_xyz == 'pdb_and_xyz':
        #use these 3 lines below and run with xyz file to speed up analsyis and not consume > 100GB memory.
        # this analysis will be an estimated 10x faster
        xyz_filename = filepath + '/' + 'SPCE_PORE_NVT_10_BOX_0_water_only.xyz'
        pdb_filename = filepath + '/' + 'SPCE_PORE_NVT_10_BOX_0_at_0M_steps_water_only.pdb'
        traj = md.load(xyz_filename, top=pdb_filename)[int(25 * 10 ** 6 / 1000):]
    elif use_pdb_only_or_pdb_and_xyz == 'pdb_only':
        # use these lines below to use the raw GOMC pdb files for analysis, but is will consume > 100GB memory.
        # this analysis will be an estimated 10x slower than using the 'pdb_and_xyz' method
        xyz_filename = filepath + '/' + 'SPCE_PORE_NVT_10_BOX_0.pdb'
        pdb_filename = filepath + '/' + 'SPCE_PORE_NVT_10_BOX_0.pdb'
        traj = md.load( xyz_filename)[int(25 * 10 ** 6 / 1000):]
    else:
        print("ERROR: neither use_pdb_only_or_pdb_and_xyz = 'pdb_only' or 'pdb_and_xyz'")

    print('traj.n_frames = '+str(traj.n_frames))

    y1, x1 = analysis.compute_s(traj,
                                surface_normal_dim=2,
                                pore_center=pore_center,
                                max_distance=0.5,
                                bin_width=0.01 )


    order_param_data = pd.DataFrame(np.column_stack([y1, x1]))

    print( order_param_data)

    order_param_data.to_csv(f'{filepath}/order_param.txt', sep="	", header=['distance_nm', 'order_param'])
    order_param_data.to_csv(f'{filepath}/order_param.csv', sep=",", header=['distance_nm', 'order_param'])
    
    plt.plot(y1, x1)
    plt.xlabel('Distance (nm)')
    plt.ylabel('S')
    plt.xlim((-0.5, 0.5))
    plt.ylim((-0.5, 0.25))
    plt.savefig(f'{filepath}/order_param.pdf')

    plt.cla()  # Clear axis
    plt.clf()  # Clear figure

Index_for_sets_list = []
distance_sets_list = []
order_param_sets_list = []

for j in range(len(filepath_list)):
    filepath = filepath_list[j]

    reading_file_Box_0 =filepath +'/'+'order_param.csv'

    Column_index_Title = ''  #
    Column_distance_Title = 'distance_nm'
    Column_order_param_Title = 'order_param'  # column title Title for iteration value


    Extracted_Data_file_Titles = [Column_index_Title, Column_distance_Title, Column_order_param_Title ]




    #*************************
    #drawing in data from single file and extracting specific rows for the liquid box (start)
    # *************************
    data_Box_0 = pd.read_csv(reading_file_Box_0, names=Extracted_Data_file_Titles, sep=',', header=0, na_values='NaN',usecols=[0,1,2], index_col=False)

    Iteration_index = data_Box_0.loc[:,Column_index_Title]
    Iteration_index = list(Iteration_index)
    Iteration_index = np.transpose(Iteration_index)
    Index_for_sets_list.append(list(Iteration_index))

    Iteration_distance = data_Box_0.loc[:, Column_distance_Title]
    Iteration_distance = list(Iteration_distance)
    Iteration_distance = np.transpose(Iteration_distance)
    distance_sets_list.append(list(Iteration_distance))

    Iteration_order_param = data_Box_0.loc[:, Column_order_param_Title]
    Iteration_order_param = list(Iteration_order_param)
    Iteration_order_param = np.transpose(Iteration_order_param)
    order_param_sets_list.append(list(Iteration_order_param))



    #*************************
    #drawing in data from single file and extracting specific rows for the vapor box (end)
    # *************************

#print(Index_for_sets_list)

Avg_Iteration_index_sets_list = []
StdDev_Iteration_index_sets_list = []

Avg_distance_sets_list = []
StdDev_distance_sets_list = []

Avg_Iteration_order_param_sets_list = []
StdDev_Iteration_order_param_sets_list = []

for h in range(0, len(Index_for_sets_list[0])):
    for g in range(0, len(Index_for_sets_list)):
        if g ==0:
            Index_for_sets_interation_list = [Index_for_sets_list[g][h]]
        else:
            Index_for_sets_interation_list.append(Index_for_sets_list[g][h])

    Avg_Iteration_index_sets_list.append(np.nanmean(Index_for_sets_interation_list))
    StdDev_Iteration_index_sets_list.append(np.std(Index_for_sets_interation_list, ddof=1))

for h in range(0, len(distance_sets_list[0])):
    for g in range(0, len(distance_sets_list)):
        if g ==0:
            distance_for_sets_interation_list = [distance_sets_list[g][h]]
        else:
            distance_for_sets_interation_list.append(distance_sets_list[g][h])

    Avg_distance_sets_list.append(np.nanmean(distance_for_sets_interation_list))
    StdDev_distance_sets_list.append(np.std(distance_for_sets_interation_list, ddof=1))

for h in range(0, len(order_param_sets_list[0])):
    for g in range(0, len(order_param_sets_list)):
        if g ==0:
            order_param_sets_interation_list = [order_param_sets_list[g][h]]
        else:
            order_param_sets_interation_list.append(order_param_sets_list[g][h])

    Avg_Iteration_order_param_sets_list.append(np.nanmean(order_param_sets_interation_list))
    StdDev_Iteration_order_param_sets_list.append(np.std(order_param_sets_interation_list, ddof=1))




Box_0_data_dataframe =pd.DataFrame(np.column_stack([Avg_distance_sets_list,
                                                    Avg_Iteration_order_param_sets_list,
                                                    StdDev_Iteration_order_param_sets_list
                                                    ]))

Box_0_data_dataframe.to_csv('Average_order_param_data.txt', sep="	",
                            header=['distance_nm', 'Avg_order_param', 'StdDev_order_param'])

Box_0_data_dataframe.to_csv('Average_order_param_data.csv', sep=",",
							header=['distance_nm', 'Avg_order_param', 'StdDev_order_param'])

plt.plot(Avg_distance_sets_list, Avg_Iteration_order_param_sets_list)
plt.xlabel('Distance (nm)')
plt.ylabel('S_order_param')
plt.xlim((-0.5, 0.5))
plt.ylim((-0.5, 0.25))
plt.savefig('Average_order_param.pdf')
