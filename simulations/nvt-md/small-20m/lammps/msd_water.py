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
def tolerant_mean(arrs):
    lens = [len(i) for i in arrs]
    arr = np.ma.empty((np.max(lens),len(arrs)))
    arr.mask = True
    for idx, l in enumerate(arrs):
        arr[:len(l),idx] = l
    return arr.mean(axis = -1), arr.std(axis=-1)


def main():
    data_path = "msd_data_water"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    os.chdir(data_path)

    top="../../../20m_small.pdb"
    project = signac.get_project()
    for (temp), group in project.groupby(("T")):
        functional="FF-UND"
        os.makedirs("{}K_{}".format(temp, functional))
        os.chdir("{}K_{}".format(temp, functional))

        traj_list=[]
        msd_list=[]
        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            discard_frames = 10000
            one_traj = md.load(job.fn("unwrap_mol.lammpstrj"), top=top)
            print("traj {} has {} ps prod".format(job, (one_traj.n_frames-discard_frames)*1))
            if one_traj.n_frames>discard_frames:
                traj_list.append(one_traj[discard_frames:])
                one_traj=one_traj[discard_frames:]
            wat_O_indices=one_traj.top.select("element O and resname wat")
            wat_xyz=one_traj.xyz[:, wat_O_indices, :]
            msd_f= freud.msd.MSD()
            msd_f.compute(wat_xyz,reset=True)
            time=msd_f.msd
            for i in range(time.shape[0]):
                time[i]=0+i*1
            msdA2=1e2*msd_f.msd
            timestep=1
            start=1000
            start_index=int(start/timestep)
            end=20000
            end_index=int(end/timestep)
            msd_fit=msdA2[start_index:end_index]
            time_fit= time[start_index:end_index]
            msd_fit=msd_fit.reshape((msd_fit.shape[0],1))
            time_fit=time_fit.reshape((msd_fit.shape[0],1))
            reg = LinearRegression().fit(time_fit, msd_fit)
            print('Coefficients: and intercept', reg.coef_, reg.intercept_)
            msd_pred=reg.predict(time_fit)
            r2_score(msd_fit, msd_pred)
            print("scikitlearn R2 score is {}".format(r2_score(msd_fit, msd_pred)))


            fig, ax = plt.subplots()

            plt.loglog(time, msdA2)
            plt.axvline(x=start,color='red',alpha=0.3, ls='--')
            plt.axvline(x=end,color='red',alpha=0.3, ls='--')
            ax.set_xlabel('Time /ps')
            ax.set_ylabel(r"MSD / $\mathrm{\AA}^2$")


            msd_list.append(msdA2)
            msd_dat  = np.column_stack([time, msdA2])
            fig.tight_layout()
            np.savetxt(
                "msd_{}.txt".format(job),
                msd_dat,
                header="time (ps)\tmsd (pm2)",
                )
            fig.savefig(
                "msd_{}.pdf".format(
                    job
                )
                )

        print(len(msd_list))
        # Now we need to do bootstrapping
        num_boots=3
        num_sample_each_boot=3
        D_values=[]
        for boot in range(num_boots):

            traj_indices=random.choices(list(range(4)), k=num_sample_each_boot)
            print("trajs chosen are {}".format(traj_indices))
            msd_boot=[]
            for index in traj_indices:
                msd_boot.append(msd_list[index])
            average_msdA2, error = tolerant_mean(msd_boot)
            time=np.copy(average_msdA2)
            for i in range(time.shape[0]):
                time[i]=0+i*1
            fig, ax = plt.subplots()

            plt.loglog(time, average_msdA2)

            plt.axvline(x=start,color='red',alpha=0.3, ls='--')
            plt.axvline(x=end,color='red',alpha=0.3, ls='--')
            ax.set_xlabel('Time /ps')
            ax.set_ylabel(r"MSD / $\mathrm{\AA}^2$")


            msd_dat  = np.column_stack([time, average_msdA2])
            fig2, ax = plt.subplots()
            log_time=np.log10(time)
            log_msdA2=np.log10(average_msdA2)
            ax.set_xlabel('log Time')
            ax.set_ylabel(r"log MSD")

            log_time_fit=log_time[start_index:end_index]
            log_msdA2_fit=log_msdA2[start_index:end_index]

            log_msdA2_fit=log_msdA2_fit.reshape((log_msdA2_fit.shape[0],1))
            log_time_fit=log_time_fit.reshape((log_time_fit.shape[0],1))
            reg = LinearRegression().fit(log_time_fit, log_msdA2_fit)
            print('Coefficients: and intercept', reg.coef_, reg.intercept_)
            log_msdA2_pred=reg.predict(log_time_fit)
            r2_score(log_msdA2_fit, log_msdA2_pred)
            print("scikitlearn R2 score is {}".format(r2_score(log_msdA2_fit, log_msdA2_pred)))
            plt.plot(log_time, log_msdA2)
            plt.plot(log_time_fit, log_msdA2_pred, label='pred')
            plt.legend()
            np.savetxt(
                "msd_boot_{}.txt".format(boot),
                msd_dat,
                header="time (ps)\tmsd (A2)",
                )
            fig.tight_layout()
            fig2.tight_layout()
            fig.savefig(
                "msd_boot_{}.pdf".format(
                    boot
                )
                )
            fig2.savefig(
                "msd_log_fit_boot_{}.pdf".format(
                    boot
                )
                )
            #Now find D
            msdA2_fit=average_msdA2[start_index:end_index]
            time_fit=time[start_index:end_index]
            msdA2_fit=msdA2_fit.reshape((msdA2_fit.shape[0],1))
            time_fit=time_fit.reshape((time_fit.shape[0],1))
            reg = LinearRegression().fit(time_fit, msdA2_fit)
            print("The diffusion coeff is {} m2/s".format(1e-20*1e12*reg.coef_/6))
            msdA2_pred=reg.predict(time_fit)
            r2_score(msdA2_fit, msdA2_pred)
            print("scikitlearn R2 score is {}".format(r2_score(msdA2_fit, msdA2_pred)))
            D_values.append(1e-20*1e12*reg.coef_/6)
            plt.close()

        D_values=np.array(D_values)
        np.savetxt(
                "diff.txt".format(boot),
                np.hstack((np.mean(D_values), np.std(D_values)/D_values.shape[0])),
                header="Diff coeff (m2/s)\tSEM",
                )
        print(D_values)
        os.chdir("..")
if __name__ == "__main__":
    main()

