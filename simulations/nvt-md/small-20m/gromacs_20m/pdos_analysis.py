import numpy as np
import signac
import shutil
from scipy import stats
import freud
import mbuild as mb
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import random
import os
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("pdf")
from glob import glob



def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def main():
    data_path = "pdos_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    top="../../10m_small.pdb"
    project = signac.get_project()
    for (temp, functional), group in project.groupby(("T","Functional")):
        print(temp,functional)
        if functional=="BLYP":
            continue
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))
        base_dir=os.getcwd()
        band_gap_list=[]
        base_dir_pdos_files= "../../pdos_rawdata/{}K_{}/".format(temp, functional)
        os.chdir(base_dir_pdos_files)
        folder_list=glob("*/")
        overall_pdos_list=[]
        water_pdos_list=[]
        TFSI_pdos_list=[]
        for folder in folder_list:
            os.chdir(folder)
            if os.path.isfile("litfsi-pdos_litfsi-k10-1.pdos"):
                pdos=np.genfromtxt("litfsi-pdos_litfsi-k10-1.pdos", skip_header=2)
                pdos[:, 2]-=1
                zero_crossings = np.where(np.diff(np.sign(pdos[:,2])))[0]
                band_gap_Ha=-(pdos[zero_crossings, 1])+(pdos[zero_crossings+1, 1])
                band_gap_ev=27.211399*band_gap_Ha
                band_gap_list.append(band_gap_ev[0])
            if os.path.isfile("overall_pdos.txt"):
                overall_pdos=np.genfromtxt("overall_pdos.txt")
                overall_pdos_list.append(overall_pdos)

            if os.path.isfile("TFSI_pdos.txt"):
                TFSI_pdos=np.genfromtxt("TFSI_pdos.txt")
                TFSI_pdos_list.append(TFSI_pdos)

            if os.path.isfile("water_pdos.txt"):
                water_pdos=np.genfromtxt("water_pdos.txt")
                water_pdos_list.append(water_pdos)


            os.chdir("..")

        os.chdir(base_dir)

        block_band_gap=list(chunks(band_gap_list,int(len(band_gap_list)/4)))
        mean=[]
        for block in block_band_gap:
            mean.append(np.mean(block))
        mean_mean=np.array([np.mean(mean)])
        mean_sem=np.array([np.std(mean)/((len(mean))**0.5)])
        mean_mean=np.reshape(mean_mean, (mean_mean.shape[0],1))
        mean_sem=np.reshape(mean_sem, (mean_sem.shape[0],1))
        print("The mean band gap is {} eV and SEM is {} eV".format(mean_mean, mean_sem))
        np.savetxt("bandgap.txt", np.hstack((mean_mean, mean_sem)), header= "band gap (eV) \t SEM")

        mean_overall_pdos=np.mean(overall_pdos_list, axis=0)
        mean_TFSI_pdos=np.mean(TFSI_pdos_list, axis=0)
        mean_water_pdos=np.mean(water_pdos_list, axis=0)
        np.savetxt("mean_overall_pdos.txt", mean_overall_pdos, header="Energy (eV) \t pdos")
        np.savetxt("mean_TFSI_pdos.txt", mean_TFSI_pdos, header="Energy (eV) \t pdos")
        np.savetxt("mean_water_pdos.txt", mean_water_pdos, header="Energy (eV) \t pdos")

        os.chdir("..")

if __name__ == "__main__":
    main()

