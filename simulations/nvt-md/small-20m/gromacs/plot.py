import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import signac
import shutil
from scipy import stats
import freud
import mbuild as mb
from scattering import scattering


def main():
    data_path = "sfac_data"
    os.chdir(data_path)
    exp_sfac_file="../298_exp_sfac.txt"
    project = signac.get_project()
    simulation_sfac={}
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    fig=plt.figure()
    exp_data_298=np.genfromtxt(exp_sfac_file,skip_header=1)
    exp_q_298=exp_data_298[:,0]
    exp_s_298=exp_data_298[:,1]
    plt.scatter(exp_q_298,exp_s_298,label="Exp.",c=u'#7f7f7f',s=8)
    for (temp, functional), group in project.groupby(("T","Functional")):
        if temp>300:
            continue
        print(temp,functional)
        os.chdir("{}K_{}".format(temp, functional))
        simulation_sfac["{}K_{}".format(temp, functional)]=np.genfromtxt('sfac.txt',skip_header=1)
        data=simulation_sfac["{}K_{}".format(temp, functional)]
        q=data[:,0]
        s=data[:,1]
        plt.plot(q,s,label="Sim. {}-D3".format( (functional)))
        plt.grid(alpha=0.3)
        plt.xlabel('$\it{q} \; \mathrm{(1/\AA)}$')
        plt.ylabel('$S(q)$')
        os.chdir("..")
    plt.xlim([0,12])
    plt.ylim([0.2,1.5])
    plt.legend()
    fig.tight_layout()
    plt.savefig('sfac_overall.png')
    plt.close()
if __name__ == "__main__":
    main()

