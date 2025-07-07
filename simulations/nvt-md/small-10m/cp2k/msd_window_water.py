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


def main():
    data_path = "msd_data_water"
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

        traj_list=[]
        msd_list=[]
        save_string = ""
        for job in group:
            seed=job.sp.Seed
            length=job.sp.L
            one_traj = md.load(job.fn("litfsi-pos-1.xyz"), top=top)
            print("traj {} has {}ps prod".format(job, (one_traj.n_frames-5000)*0.005))
            if one_traj.n_frames>5000:
                traj_list.append(one_traj[5000:])
                one_traj=one_traj[5000:]
            wat_O_indices=one_traj.top.select("element O and resname wat")
            wat_xyz=one_traj.xyz[:, wat_O_indices, :]

            window_width = 70
            skip_time= 1
            timestep=0.005
            list_msds=[]
            msd_window = 40  #finding a window of slope 1


            skip_steps = int(skip_time/timestep)
            window_width_steps = int(window_width/timestep)
            i=0
            while (i*skip_steps+window_width_steps) < len(wat_xyz):
                print(i*skip_steps+window_width_steps)
                window_xyz = wat_xyz[i*skip_steps:(i*skip_steps+window_width_steps)]
                coords0 = window_xyz[0]
                window_msdA2 = []

                for frame in range(len(window_xyz)):
                    diff = window_xyz[frame]-coords0
                    window_msdA2.append(100*np.mean(np.linalg.norm(diff, axis=1)**2))
                window_msdA2= np.array(window_msdA2)
                list_msds.append(window_msdA2)
                i+=1
            msdA2= np.mean(list_msds, axis=0)
            timestep = 0.005 # this needs to be the actual time between frames
            time = np.arange(len(msdA2))*timestep # make the lag-time axis





            fig, ax = plt.subplots()

            plt.loglog(time, msdA2)
            #plt.axvline(x=start,color='red',alpha=0.3, ls='--')
            #plt.axvline(x=end,color='red',alpha=0.3, ls='--')
            ax.set_xlabel('Time /ps')
            ax.set_ylabel('MSD / A$^2$')
            plt.tight_layout()
            msd_list.append(msdA2)
            msd_dat  = np.column_stack([time, msdA2])

            np.savetxt(
                "msd_{}.txt".format(job),
                msd_dat,
                header="time (ps)\tmsd (pm2)",
                )
            #fig.savefig("msd_{}.pdf".format( job ) )

        print("We have {} samples for msd".format(len(msd_list)))
        # Now we need to do bootstrapping
        num_boots=1
        num_sample_each_boot=4
        D_values=[]
        for boot in range(num_boots):

            #traj_indices=random.choices(list(range(4)), k=num_sample_each_boot)
            traj_indices = [0,1,2,3]
            print("trajs chosen are {}".format(traj_indices))
            msd_boot=[]
            for index in traj_indices:
                msd_boot.append(msd_list[index])
            average_msdA2 = np.mean(msd_boot, axis =0)
            print(average_msdA2)
            time=np.copy(average_msdA2)
            for i in range(time.shape[0]):
                time[i]=0+i*0.005

            msd_dat  = np.column_stack([time, average_msdA2])
            log_time=np.log10(time)
            log_msdA2=np.log10(average_msdA2)



            orig_start_time = 10
            step_start_time = 1
            start_time = 0
            slopes = []
            start_times = []
            i = 0
            while  (start_time+msd_window ) < time[-1]-10:
                start_time = orig_start_time+i*step_start_time
                end_time = start_time+msd_window
                start_index = int(start_time/timestep)
                end_index = int(end_time/timestep)
                print("start and end indices are {}, {}".format(start_index, end_index))

                log_time_fit       =    log_time[start_index:end_index]
                log_msdA2_fit      =    log_msdA2[start_index:end_index]

                log_msdA2_fit=log_msdA2_fit.reshape((log_msdA2_fit.shape[0],1))
                log_time_fit=log_time_fit.reshape((log_time_fit.shape[0],1))


                reg = LinearRegression().fit(log_time_fit, log_msdA2_fit)
                print("Log log slope when start_time = {} ps and msd_window = {} ps is {}".format(start_time, msd_window, reg.coef_))
                slopes.append(reg.coef_)
                start_times.append(start_time)
                i=i+1

            abs_slope_diff = [abs(k-1) for k in slopes]
            index_min = min(range(len(abs_slope_diff)), key=abs_slope_diff.__getitem__) 
    
            print("Minimum slope is {} and has a start_time value of {}".format(slopes[index_min], start_times[index_min]))
            print("Slopes are {}\n Start times are {}\n".format(slopes, start_times))
            save_string+="Minimum slope is {} and has a start_time value of {}\n".format(slopes[index_min], start_times[index_min])
            save_string+="Slopes are {}\n Start times are {}\n".format(slopes, start_times)

            start_time = start_times[index_min]
            end_time   = start_times[index_min] + msd_window

            start_index = int(start_time/timestep)
            end_index = int(end_time/timestep)

            fig, ax = plt.subplots()

            plt.loglog(time, average_msdA2)

            ax.set_xlabel('Time /ps')
            ax.set_ylabel('MSD / A$^2$')
            plt.axvline(x=start_time,color='red',alpha=0.3, ls='--')
            plt.axvline(x=end_time,color='red',alpha=0.3, ls='--')
            plt.tight_layout()

            log_time_fit=log_time[start_index:end_index]
            log_msdA2_fit=log_msdA2[start_index:end_index]

            log_msdA2_fit=log_msdA2_fit.reshape((log_msdA2_fit.shape[0],1))
            log_time_fit=log_time_fit.reshape((log_time_fit.shape[0],1))
            reg = LinearRegression().fit(log_time_fit, log_msdA2_fit)
            print('Coefficients: and intercept', reg.coef_, reg.intercept_)
            log_msdA2_pred=reg.predict(log_time_fit)
            r2_score(log_msdA2_fit, log_msdA2_pred)
            print("scikitlearn R2 score is {}".format(r2_score(log_msdA2_fit, log_msdA2_pred)))
            fig2, ax = plt.subplots()
            plt.plot(log_time, log_msdA2)
            plt.plot(log_time_fit, log_msdA2_pred, label='pred')
            plt.legend()
            np.savetxt(
                "msd_boot_{}.txt".format(boot),
                msd_dat,
                header="time (ps)\tmsd (A2)",
                )
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
            plt.close()

            fig3, ax = plt.subplots()
            plt.plot(time, average_msdA2)

            ax.set_xlabel('Time /ps')
            ax.set_ylabel('MSD / A$^2$')
            plt.axvline(x=start_time,color='red',alpha=0.3, ls='--')
            plt.axvline(x=end_time,color='red',alpha=0.3, ls='--')
            plt.tight_layout()
            fig3.savefig(
                "msd_linear_boot_{}.pdf".format(
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
        with open("diff_information.txt", "w") as text_file:
            text_file.write(save_string)

        print(D_values)
        os.chdir("..")
if __name__ == "__main__":
    main()

