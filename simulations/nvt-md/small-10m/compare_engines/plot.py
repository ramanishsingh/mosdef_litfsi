import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("pdf")
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import StrMethodFormatter, NullFormatter, MultipleLocator
import numpy as np
import string

##########################################################################################################
# Plot 1 Comparing sfacs of cp2k, lammps
##########################################################################################################
fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6,6))
ax1,ax2,ax3,ax4=axs[0,0],axs[0,1],axs[1,0],axs[1,1]
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.suptitle('Structure factor comparisoin between CP2K, LAMMPS, and GROMACS')
sim_298K_10m_cp2k_PBE=np.genfromtxt('../cp2k/sfac_data/298.0K_PBE/sfac.txt', skip_header=1)
sim_298K_10m_lam     = np.genfromtxt('../lammps/sfac_data/298.0K_FF-UND/sfac.txt', skip_header=1) 


ax1.plot(sim_298K_10m_cp2k_PBE[:,0], sim_298K_10m_cp2k_PBE[:,1], label='Sim PBE')
ax1.plot(sim_298K_10m_lam [:,0], sim_298K_10m_lam [:,1], label='Exp', linestyle = 'None', marker='s', markersize=0.5)

#adding text
ax1.text(0.83, 0.3, '298 K 10m', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
ax2.text(0.83, 0.3, '298 K 20m', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
ax3.text(0.83, 0.3, '373 K 10m', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes)
ax4.text(0.83, 0.3, '373 K 20m', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes)
ax1.legend(frameon=True, fontsize=9, labelspacing=0.05)
ax3.set_xlim([0,8])
ax4.set_xlim([0,8])
ax1.set_ylim([0,1.5])
ax3.set_ylim([0,1.5])
ax3.set_xlabel(r"$q$/$\mathrm{\AA}^{-1}$")
ax4.set_xlabel(r"$q$/$\mathrm{\AA}^{-1}$")
ax1.set_ylabel("$S(q)$")
ax3.set_ylabel("$S(q)$")
#a b c d
axs = axs.flat
for n, ax in enumerate(axs):
    ax.text(-0.0, 1.05, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,)
fig.tight_layout()
fig.savefig("sfac.pdf")
plt.close(fig)
'''
##########################################################################################################
##########################################################################################################


##########################################################################################################
# Plot 2 , individual sfac for 298 K 10 m and 20 m
########################################################################################################## 
fig, axs = plt.subplots(1, 2,sharey=True, figsize=(6,3.5))

ax1,ax2=axs[0],axs[1]
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.suptitle('Structure factor decomposition')

LiLi_298K_10m      =    np.genfromtxt('../small-10m/cp2k/sfac_data/298.0K_PBE/sfac_Li-Li.txt', skip_header=1)
LiTFSI_298K_10m    =    np.genfromtxt('../small-10m/cp2k/sfac_data/298.0K_PBE/sfac_Li-TFSI.txt', skip_header=1)
LiH2O_298K_10m     =    np.genfromtxt('../small-10m/cp2k/sfac_data/298.0K_PBE/sfac_Li-H2O.txt', skip_header=1)
TFSITFSI_298K_10m  =    np.genfromtxt('../small-10m/cp2k/sfac_data/298.0K_PBE/sfac_TFSI-TFSI.txt', skip_header=1)
TFSIH2O_298K_10m   =    np.genfromtxt('../small-10m/cp2k/sfac_data/298.0K_PBE/sfac_TFSI-H2O.txt', skip_header=1)
H2OH2O_298K_10m    =    np.genfromtxt('../small-10m/cp2k/sfac_data/298.0K_PBE/sfac_H2O-H2O.txt', skip_header=1)


LiLi_298K_20m      =    np.genfromtxt('../small-20m/cp2k/sfac_data/298.0K_PBE/sfac_Li-Li.txt', skip_header=1)
LiTFSI_298K_20m    =    np.genfromtxt('../small-20m/cp2k/sfac_data/298.0K_PBE/sfac_Li-TFSI.txt', skip_header=1)
LiH2O_298K_20m     =    np.genfromtxt('../small-20m/cp2k/sfac_data/298.0K_PBE/sfac_Li-H2O.txt', skip_header=1)
TFSITFSI_298K_20m  =    np.genfromtxt('../small-20m/cp2k/sfac_data/298.0K_PBE/sfac_TFSI-TFSI.txt', skip_header=1)
TFSIH2O_298K_20m   =    np.genfromtxt('../small-20m/cp2k/sfac_data/298.0K_PBE/sfac_TFSI-H2O.txt', skip_header=1)
H2OH2O_298K_20m    =    np.genfromtxt('../small-20m/cp2k/sfac_data/298.0K_PBE/sfac_H2O-H2O.txt', skip_header=1)

ax1.plot(LiLi_298K_10m[:,0], LiLi_298K_10m[:,1], label = 'Li-Li')
ax1.plot(LiTFSI_298K_10m[:,0], LiTFSI_298K_10m[:,1], label='Li-TFSI')
ax1.plot(LiH2O_298K_10m[:,0], LiH2O_298K_10m[:,1], label='Li-H$_2$O')
ax1.plot(TFSITFSI_298K_10m[:,0], TFSITFSI_298K_10m[:,1], label='TFSI-TFSI')
ax1.plot(TFSIH2O_298K_10m[:,0], TFSIH2O_298K_10m[:,1], label='TFSI-H$_2$O')
ax1.plot(H2OH2O_298K_10m[:,0], H2OH2O_298K_10m[:,1], label= 'H$_2$O-H$_2$O')

ax2.plot(LiLi_298K_20m[:,0], LiLi_298K_20m[:,1], label = 'Li-Li')
ax2.plot(LiTFSI_298K_20m[:,0], LiTFSI_298K_20m[:,1], label='Li-TFSI')
ax2.plot(LiH2O_298K_20m[:,0], LiH2O_298K_20m[:,1], label='Li-H$_2$O')
ax2.plot(TFSITFSI_298K_20m[:,0], TFSITFSI_298K_20m[:,1], label='TFSI-TFSI')
ax2.plot(TFSIH2O_298K_20m[:,0], TFSIH2O_298K_20m[:,1], label='TFSI-H$_2$O')
ax2.plot(H2OH2O_298K_20m[:,0], H2OH2O_298K_20m[:,1], label= 'H$_2$O-H$_2$O')

#adding text
ax1.text(0.83, 0.3, '298 K 10m', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
ax2.text(0.83, 0.3, '298 K 20m', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
ax1.legend(frameon=False, loc='best', ncol=2, fontsize=9, labelspacing=0.05)
ax1.set_xlim([0.5,3])
ax2.set_xlim([0.5,3])
ax1.set_ylim([-3.5,1.5])
#ax3.set_ylim([0,1.5])
ax1.set_xlabel(r"$q$/$\mathrm{\AA}^{-1}$")
ax2.set_xlabel(r"$q$/$\mathrm{\AA}^{-1}$")
ax1.set_ylabel("$S(q)$")
#a b c d
axs = axs.flat
for n, ax in enumerate(axs):
    ax.text(-0.0, 1.05, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,)
fig.tight_layout()
fig.savefig("sfac_decomposed.pdf")
plt.close(fig)
##########################################################################################################
##########################################################################################################
'''

##########################################################################################################
# Plot 3 , All TFSI rdfs
##########################################################################################################
fig, axs = plt.subplots(2, 3, sharex=True, sharey='row', figsize=(9,6))

ax1,ax2,ax3,ax4,ax5,ax6=axs[0,0],axs[1,0],axs[0,1], axs[1,1],axs[0,2],axs[1,2]
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.suptitle('TFSI RDF')

TFSI_TFSI_298K_10m_cp2k_pbe      =    np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/TFSI-TFSI.txt', skip_header=1)
CTFSI_CTFSI_298K_10m_cp2k_pbe    =    np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/C(TFSI)-C(TFSI).txt', skip_header=1)
TFSI_water_298K_10m_cp2k_pbe     =    np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/TFSI-water.txt', skip_header=1)

TFSI_TFSI_298K_10m_lam      =    np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/TFSI-TFSI.txt', skip_header=1)
CTFSI_CTFSI_298K_10m_lam    =    np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/C(TFSI)-C(TFSI).txt', skip_header=1)
TFSI_water_298K_10m_lam     =    np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/TFSI-water.txt', skip_header=1)


ax1.plot(TFSI_TFSI_298K_10m_cp2k_pbe[:,0], TFSI_TFSI_298K_10m_cp2k_pbe[:,1], label='298 K 10m')
ax1.plot(TFSI_TFSI_298K_10m_lam[:,0], TFSI_TFSI_298K_10m_lam[:,1], label='298 K 20m')

ax2.plot(TFSI_TFSI_298K_10m_cp2k_pbe[:,2], TFSI_TFSI_298K_10m_cp2k_pbe[:,3], label='298 K 10m')
ax2.plot(TFSI_TFSI_298K_10m_lam[:,2], TFSI_TFSI_298K_10m_lam[:,3], label='298 K 20m')

ax3.plot(CTFSI_CTFSI_298K_10m_cp2k_pbe[:,0], CTFSI_CTFSI_298K_10m_cp2k_pbe[:,1], label='298 K 10m')
ax3.plot(CTFSI_CTFSI_298K_10m_lam[:,0], CTFSI_CTFSI_298K_10m_lam[:,1], label='298 K 20m')

ax4.plot(CTFSI_CTFSI_298K_10m_cp2k_pbe[:,2], CTFSI_CTFSI_298K_10m_cp2k_pbe[:,3], label='298 K 10m')
ax4.plot(CTFSI_CTFSI_298K_10m_lam[:,2] , CTFSI_CTFSI_298K_10m_lam[:,3], label='298 K 20m')

ax5.plot(TFSI_water_298K_10m_cp2k_pbe[:,0], TFSI_water_298K_10m_cp2k_pbe[:,1], label='298 K 10m')
ax5.plot(TFSI_water_298K_10m_lam[:,0], TFSI_water_298K_10m_lam[:,1], label='298 K 20m')

ax6.plot(TFSI_water_298K_10m_cp2k_pbe[:,2], TFSI_water_298K_10m_cp2k_pbe[:,3], label='298 K 10m')
ax6.plot(TFSI_water_298K_10m_lam[:,2], TFSI_water_298K_10m_lam[:,3], label='298 K 20m')
ax1.legend(frameon=False, loc='best', ncol=2, fontsize=9, labelspacing=0.05)
ax2.set_xlabel(r"$r$/$\mathrm{\AA}$")
ax4.set_xlabel(r"$r$/$\mathrm{\AA}$")
ax6.set_xlabel(r"$r$/$\mathrm{\AA}$")
ax1.set_ylabel(r"$g(r)$")
ax2.set_ylabel(r"$g(r)$")
ax1.set_ylim([0,2.5])
ax1.set_xlim([3,8])
ax2.set_ylim([0,10])
axs = axs.flat
for n, ax in enumerate(axs):
    ax.text(-0.0, 1.05, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,)

#adding text
ax1.text(0.8, 0.1, 'TFSI-TFSI', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
ax2.text(0.2, 0.1, 'TFSI-TFSI', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
ax3.text(0.8, 0.1, 'C(TFSI)-C(TFSI)', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes)
ax4.text(0.8, 0.1, 'C(TFSI)-C(TFSI)', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes)
ax5.text(0.8, 0.1, 'TFSI-water', horizontalalignment='center', verticalalignment='center', transform=ax5.transAxes)
ax6.text(0.8, 0.1, 'TFSI-water', horizontalalignment='center', verticalalignment='center', transform=ax6.transAxes)


fig.tight_layout()


fig.savefig("TFSI_RDF.pdf")
plt.close(fig)
##########################################################################################################
##########################################################################################################

##########################################################################################################
# Plot 4 , All water rdfs and distribution
##########################################################################################################
fig, axs = plt.subplots(2, 2, figsize=(6,6))

ax1,ax2,ax3,ax4=axs[0,0],axs[0,1],axs[1,0], axs[1,1]
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.suptitle('Water RDF and Cluster size distribution')

Owater_Owater_298K_10m_cp2k_PBE  =    np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/O(water)-O(water).txt', skip_header=1)
Owater_Owater_298K_10m_lam  =    np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/O(water)-O(water).txt', skip_header=1)

prob_Owater_Owater_298K_10m_cp2k_PBE  = np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/prob_O(water)_O(water)_cluster.txt', skip_header=1)
prob_Owater_Owater_298K_10m_lam  = np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/prob_O(water)_O(water)_cluster.txt', skip_header=1)

cum_prob_O_O_298K_10m_cp2k_PBE  =  np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/cum_prob_O_O.txt', skip_header=1)
cum_prob_O_O_298K_10m_lam  =  np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/cum_prob_O_O.txt', skip_header=1)

ax1.plot(Owater_Owater_298K_10m_cp2k_PBE[:,0], Owater_Owater_298K_10m_cp2k_PBE[:,1], label='298 K 10m')
ax1.plot(Owater_Owater_298K_10m_lam[:,0], Owater_Owater_298K_10m_lam[:,1], label='298 K 20m')

ax2.plot(Owater_Owater_298K_10m_cp2k_PBE[:,2], Owater_Owater_298K_10m_cp2k_PBE[:,3], label='298 K 10m')
ax2.plot(Owater_Owater_298K_10m_lam[:,2], Owater_Owater_298K_10m_lam[:,3], label='298 K 20m')

ax3.plot(prob_Owater_Owater_298K_10m_cp2k_PBE[:,0], prob_Owater_Owater_298K_10m_cp2k_PBE[:,1], '-o',label='298 K 10m')
ax3.plot(prob_Owater_Owater_298K_10m_lam[:,0], prob_Owater_Owater_298K_10m_lam[:,1], '-s',label='298 K 20m')

ax4.loglog(cum_prob_O_O_298K_10m_cp2k_PBE[:,0], cum_prob_O_O_298K_10m_cp2k_PBE[:,1], label='298 K 10m')
ax4.loglog(cum_prob_O_O_298K_10m_lam[:,0], cum_prob_O_O_298K_10m_lam[:,1], label='298 K 20m')

ax1.legend(frameon=False, loc='best', ncol=1, fontsize=9, labelspacing=0.05)
ax1.set_xlabel(r"$r$/$\mathrm{\AA}$")
ax2.set_xlabel(r"$r$/$\mathrm{\AA}$")
ax3.set_xlabel("O-O Coordination number")
ax4.set_xlabel("Fractional cluster size")
ax1.set_ylabel(r"$g(r)$")
ax2.set_ylabel(r"$n(r)$")
ax3.set_ylabel("Probability")
ax4.set_ylabel("Cumulative probability")

ax1.set_xlim([2,8])
ax2.set_xlim([2,8])
ax3.set_xlim([-0.5,6.5])
ax3.xaxis.set_major_locator(MaxNLocator(integer=True))
ax3.minorticks_off()
ax4.set_xlim([1e-2,1])
ax1.set_ylim([0,4.5])
ax2.set_ylim([0,10])
ax3.set_ylim([0,0.4])
ax4.set_ylim([1e-1,1])
ax4.yaxis.set_major_locator(MultipleLocator(0.5))
ax4.yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax4.yaxis.set_minor_formatter(NullFormatter())
#ax4.xaxis.set_major_locator(MultipleLocator(10))
ax4.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax4.xaxis.set_minor_formatter(NullFormatter())
axs = axs.flat
for n, ax in enumerate(axs):
    ax.text(-0.0, 1.05, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,)

#adding text
ax1.text(0.73, 0.1, 'O(water)-O(water)', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
ax2.text(0.73, 0.1, 'O(water)-O(water)', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
ax3.text(0.73, 0.9, 'O(water)-O(water)', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes)
ax4.text(0.73, 0.1, 'O(water)-O(water)', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes)


fig.tight_layout()


fig.savefig("water_RDF.pdf")
plt.close(fig)
##########################################################################################################
##########################################################################################################

##########################################################################################################
# Plot 5 , Li RDF
##########################################################################################################
fig, axs = plt.subplots(3, 2, figsize=(6,9))

ax1,ax2,ax3,ax4,ax5,ax6=axs[0,0],axs[0,1],axs[1,0], axs[1,1],axs[2,0],axs[2,1]
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.suptitle('Li RDF and coordination')

Li_OTFSI_298K_10m_cp2k_PBE       =    np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/Li-O(TFSI).txt', skip_header=1)
Li_OTFSI_298K_10m_lam       =    np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/Li-O(TFSI).txt', skip_header=1)


Li_Owater_298K_10m_cp2k_PBE      =    np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/Li-O(water).txt', skip_header=1)
Li_Owater_298K_10m_lam      =    np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/Li-O(water).txt', skip_header=1)

prob_Li_Owater_298K_10m_cp2k_PBE =    np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/prob_Li_O(water)_cluster.txt', skip_header=1)
prob_Li_Owater_298K_10m_lam =    np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/prob_Li_O(water)_cluster.txt', skip_header=1)

prob_Li_OTFSI_298K_10m_cp2k_PBE  =    np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/prob_Li_O(TF2)_cluster.txt', skip_header=1)
prob_Li_OTFSI_298K_10m_lam  =    np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/prob_Li_O(TF2)_cluster.txt', skip_header=1)


ax1.plot(Li_OTFSI_298K_10m_cp2k_PBE[:,0], Li_OTFSI_298K_10m_cp2k_PBE[:,1], label='298 K 10m')
ax1.plot(Li_OTFSI_298K_10m_lam[:,0], Li_OTFSI_298K_10m_lam[:,1], label='298 K 20m')


ax3.plot(Li_OTFSI_298K_10m_cp2k_PBE[:,2], Li_OTFSI_298K_10m_cp2k_PBE[:,3], label='298 K 10m')
ax3.plot(Li_OTFSI_298K_10m_lam[:,2], Li_OTFSI_298K_10m_lam[:,3], label='298 K 20m')


ax2.plot(Li_Owater_298K_10m_cp2k_PBE[:,0], Li_Owater_298K_10m_cp2k_PBE[:,1], label='298 K 10m')
ax2.plot(Li_Owater_298K_10m_lam[:,0], Li_Owater_298K_10m_lam[:,1], label='298 K 20m')

ax4.plot(Li_Owater_298K_10m_cp2k_PBE[:,2], Li_Owater_298K_10m_cp2k_PBE[:,3], label='298 K 10m')
ax4.plot(Li_Owater_298K_10m_lam[:,2], Li_Owater_298K_10m_lam[:,3], label='298 K 20m')



ax5.plot(prob_Li_OTFSI_298K_10m_cp2k_PBE[:,0], prob_Li_OTFSI_298K_10m_cp2k_PBE[:,1], '-o', label='298 K 10m')
ax5.plot(prob_Li_OTFSI_298K_10m_lam[:,0], prob_Li_OTFSI_298K_10m_lam[:,1], '-s', label='298 K 20m')


ax6.plot(prob_Li_Owater_298K_10m_cp2k_PBE[:,0], prob_Li_Owater_298K_10m_cp2k_PBE[:,1],  '-o',label='298 K 10m')
ax6.plot(prob_Li_Owater_298K_10m_lam[:,0], prob_Li_Owater_298K_10m_lam[:,1],  '-s', label='298 K 20m')

ax1.legend(frameon=False, loc='best', ncol=2, fontsize=7, labelspacing=0.05)

#axes lables
ax1.get_shared_y_axes().join(ax1, ax2)
ax1.get_shared_x_axes().join(ax1, ax3)
ax3.get_shared_y_axes().join(ax3, ax4)
ax2.get_shared_x_axes().join(ax2, ax4)
ax5.get_shared_y_axes().join(ax5, ax6)
ax5.get_shared_x_axes().join(ax5, ax6)

ax3.set_xlabel(r"$r$/$\mathrm{\AA}$")
ax4.set_xlabel(r"$r$/$\mathrm{\AA}$")
ax5.set_xlabel(r"Number of oxygen atoms")
ax6.set_xlabel(r"Number of oxygen atoms")
ax1.set_ylabel(r"$g(r)$")
ax3.set_ylabel(r"$n(r)$")
ax5.set_ylabel(r"Probability")

#axes limits
ax1.set_ylim([0,5])

ax1.set_xlim([1,8])
ax2.set_xlim([1,8])
ax3.set_ylim([0,10])
ax5.set_ylim([0,0.6])
ax5.set_xlim([-0.5, 6.5])

#Only integer ticks in prob
ax5.xaxis.set_major_locator(MaxNLocator(integer=True))
ax5.tick_params(axis='x', which='minor', length=0)
ax6.xaxis.set_major_locator(MaxNLocator(integer=True))
ax6.tick_params(axis='x', which='minor', length=0)

ax1.legend(frameon=False, loc='best', fontsize=9, labelspacing=0.05)

axs = axs.flat
for n, ax in enumerate(axs):
    ax.text(-0.0, 1.05, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,)

#adding text
ax1.text(0.83, 0.6, 'Li-O(TFSI)', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
ax2.text(0.83, 0.6, 'Li-O(water)', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
ax3.text(0.83, 0.6, 'Li-O(TFSI)', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes)
ax4.text(0.83, 0.6, 'Li-O(water)', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes)
ax5.text(0.83, 0.6, 'Li-O(TFSI)', horizontalalignment='center', verticalalignment='center', transform=ax5.transAxes)
ax6.text(0.83, 0.6, 'Li-O(water)', horizontalalignment='center', verticalalignment='center', transform=ax6.transAxes)


fig.tight_layout()


fig.savefig("Li_RDF.pdf")
plt.close(fig)
##########################################################################################################
##########################################################################################################

##########################################################################################################
# Plot 6 , TFSI cumulative prob
##########################################################################################################
fig, axs = plt.subplots(1, 1, figsize=(3,3))

ax1=axs
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.suptitle('TFSI fractional cluster size')


cum_prob_S_S_298K_10m_cp2k_PBE  =  np.genfromtxt('../cp2k/rdf_data/298.0K_PBE/cum_prob_S_S.txt', skip_header=1)
cum_prob_S_S_298K_10m_lam  =  np.genfromtxt('../lammps/rdf_data/298.0K_FF-UND/cum_prob_S_S.txt', skip_header=1)
ax1.loglog(cum_prob_S_S_298K_10m_cp2k_PBE[:,0], cum_prob_S_S_298K_10m_cp2k_PBE[:,1], label='298 K 10m')
ax1.loglog(cum_prob_S_S_298K_10m_lam[:,0], cum_prob_S_S_298K_10m_lam[:,1], label='298 K 20m')
ax1.legend(frameon=False, loc='best', ncol=2, fontsize=7, labelspacing=0.05)
ax1.set_xlabel("Fractional cluster size")
ax1.set_ylabel("Cumulative probability")

ax1.set_xlim([1e-1,1])
ax1.set_ylim([2e-1,1])
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax1.yaxis.set_minor_formatter(NullFormatter())
#ax4.xaxis.set_major_locator(MultipleLocator(10))
ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax1.xaxis.set_minor_formatter(NullFormatter())
fig.tight_layout()


fig.savefig("TFSI_cumulative_prob.pdf")
plt.close(fig)
##########################################################################################################
##########################################################################################################

'''
##########################################################################################################
# Plot 7 , bandgap
##########################################################################################################
fig, axs = plt.subplots(1, 1, figsize=(3,3))

ax1=axs
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.suptitle('Band gap')

bandgap_298K_10m       =  np.genfromtxt('../small-10m/cp2k/pdos_data/298.0K_PBE/bandgap.txt', skip_header=1)
bandgap_298K_20m       =  np.genfromtxt('../small-20m/cp2k/pdos_data/298.0K_PBE/bandgap.txt', skip_header=1)
bandgap_373K_10m       =  np.genfromtxt('../small-10m/cp2k/pdos_data/373.0K_PBE/bandgap.txt', skip_header=1)
bandgap_373K_20m       =  np.genfromtxt('../small-20m/cp2k/pdos_data/373.0K_PBE/bandgap.txt', skip_header=1)

ax1.errorbar([10,20], [bandgap_298K_10m[0],bandgap_298K_20m[0]], yerr=[bandgap_298K_10m[1],bandgap_298K_20m[1]], label="298 K", fmt='-o')
ax1.errorbar([10,20], [bandgap_373K_10m[0],bandgap_373K_20m[0]], yerr=[bandgap_373K_10m[1],bandgap_373K_20m[1]], label="298 K", fmt='-o')

ax1.legend(frameon=False, loc='best', ncol=2, fontsize=7, labelspacing=0.05)
ax1.set_xlabel("LiTFSI molality")
ax1.set_ylabel("Band gap (eV)")

#ax1.set_xlim([1e-1,1])
#ax1.set_ylim([2e-1,1])
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax1.yaxis.set_minor_formatter(NullFormatter())
#ax4.xaxis.set_major_locator(MultipleLocator(10))
ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax1.xaxis.set_minor_formatter(NullFormatter())
fig.tight_layout()


fig.savefig("band_gap.pdf")
plt.close(fig)
##########################################################################################################
##########################################################################################################

##########################################################################################################
# Plot 5 , PDOS
##########################################################################################################
fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6,6))

ax1,ax2,ax3,ax4=axs[0,0],axs[0,1],axs[1,0], axs[1,1]
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.suptitle('Projected density of states')

water_298K_10m          =    np.genfromtxt('../small-10m/cp2k/pdos_data/298.0K_PBE/mean_water_pdos.txt', skip_header=1)
TFSI_298K_10m           =    np.genfromtxt('../small-10m/cp2k/pdos_data/298.0K_PBE/mean_TFSI_pdos.txt', skip_header=1)
overall_298K_10m        =    np.genfromtxt('../small-10m/cp2k/pdos_data/298.0K_PBE/mean_overall_pdos.txt', skip_header=1)

water_298K_20m          =    np.genfromtxt('../small-20m/cp2k/pdos_data/298.0K_PBE/mean_water_pdos.txt', skip_header=1)
TFSI_298K_20m           =    np.genfromtxt('../small-20m/cp2k/pdos_data/298.0K_PBE/mean_TFSI_pdos.txt', skip_header=1)
overall_298K_20m        =    np.genfromtxt('../small-20m/cp2k/pdos_data/298.0K_PBE/mean_overall_pdos.txt', skip_header=1)

water_373K_10m          =    np.genfromtxt('../small-10m/cp2k/pdos_data/298.0K_PBE/mean_water_pdos.txt', skip_header=1)
TFSI_373K_10m           =    np.genfromtxt('../small-10m/cp2k/pdos_data/298.0K_PBE/mean_TFSI_pdos.txt', skip_header=1)
overall_373K_10m        =    np.genfromtxt('../small-10m/cp2k/pdos_data/298.0K_PBE/mean_overall_pdos.txt', skip_header=1)

water_373K_20m          =    np.genfromtxt('../small-20m/cp2k/pdos_data/373.0K_PBE/mean_water_pdos.txt', skip_header=1)
TFSI_373K_20m           =    np.genfromtxt('../small-20m/cp2k/pdos_data/373.0K_PBE/mean_TFSI_pdos.txt', skip_header=1)
overall_373K_20m        =    np.genfromtxt('../small-20m/cp2k/pdos_data/373.0K_PBE/mean_overall_pdos.txt', skip_header=1)


ax1.plot(  water_298K_10m[:,0],  water_298K_10m[:,1], label='H$_2$O')
ax1.plot(   TFSI_298K_10m[:,0],   TFSI_298K_10m[:,1], label='TFSI')
ax1.plot(overall_298K_10m[:,0],overall_298K_10m[:,1], label='Overall')

ax2.plot(  water_373K_10m[:,0],  water_373K_10m[:,1], label='H$_2$O')
ax2.plot(   TFSI_373K_10m[:,0],   TFSI_373K_10m[:,1], label='TFSI')
ax2.plot(overall_373K_10m[:,0],overall_373K_10m[:,1], label='Overall')

ax3.plot(  water_298K_20m[:,0],  water_298K_20m[:,1], label='H$_2$O')
ax3.plot(   TFSI_298K_20m[:,0],   TFSI_298K_20m[:,1], label='TFSI')
ax3.plot(overall_298K_20m[:,0],overall_298K_20m[:,1], label='Overall')

ax4.plot(  water_373K_20m[:,0],  water_373K_20m[:,1], label='H$_2$O')
ax4.plot(   TFSI_373K_20m[:,0],   TFSI_373K_20m[:,1], label='TFSI')
ax4.plot(overall_373K_20m[:,0],overall_373K_20m[:,1], label='Overall')



ax1.legend(frameon=False, loc='best', ncol=1, fontsize=9, labelspacing=0.05)
ax3.set_xlabel("Energy / eV")
ax4.set_xlabel("Energy / eV")
ax1.set_ylabel("PDOS")
ax3.set_ylabel("PDOS")

ax1.set_xlim([-1,10])
ax2.set_xlim([-1,8])
#ax3.set_xlim([-0.5,6.5])
#ax3.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax3.minorticks_off()
#ax4.set_xlim([1e-2,1])
#ax1.set_ylim([0,4.5])
#ax2.set_ylim([0,10])
#ax3.set_ylim([0,0.4])
#ax4.set_ylim([1e-1,1])
#ax4.yaxis.set_major_locator(MultipleLocator(0.5))
#ax4.yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
#ax4.yaxis.set_minor_formatter(NullFormatter())
#ax4.xaxis.set_major_locator(MultipleLocator(10))
#ax4.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
#ax4.xaxis.set_minor_formatter(NullFormatter())
axs = axs.flat
for n, ax in enumerate(axs):
    ax.text(-0.0, 1.05, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,)

#adding text
ax1.text(0.73, 0.1, '298 K, 10m', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
ax2.text(0.73, 0.1, '373 K, 10m', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
ax3.text(0.73, 0.9, '298 K, 20m', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes)
ax4.text(0.73, 0.1, '373 K, 20m', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes)


fig.tight_layout()


fig.savefig("pdos.pdf")
plt.close(fig)
##########################################################################################################
##########################################################################################################
'''

