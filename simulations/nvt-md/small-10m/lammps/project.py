import flow
from flow import FlowProject, directives
import warnings
import shutil
import os
from subprocess import Popen, PIPE

warnings.filterwarnings("ignore", category=DeprecationWarning)


class Project(FlowProject):
    pass


@Project.label
def has_data_file(job):
    return job.isfile("data.mol")


@Project.label
def has_control_file(job):
    return job.isfile("control.inp")


@Project.label
def md_completed(job):
    "Verify that the md simulation has completed"

    try:
        return job.isfile(job.doc.output_filename)
    except:
        return False



@Project.operation
@Project.post(has_data_file)
def copy_data_file(job):
    with job:
        seed=int(job.sp.Seed)
        T=int(job.sp.T)
        shutil.copyfile(Project().root_directory() + "/lammps_input_file/"+"data_10m_small_{}K_{}".format(T,seed), job.workspace() + "/data.mol")

@Project.operation
@Project.pre(has_data_file)
@Project.post(has_control_file)
def copy_control_file(job):
    with job:
        seed=int(job.sp.Seed)
        T=int(job.sp.T)
        shutil.copyfile(Project().root_directory() + "/lammps_input_file/"+"in.{}".format(T), job.workspace() + "/control.inp")

@Project.operation
@Project.pre(has_control_file)
@Project.pre(has_data_file)
@Project.post(md_completed)
def run_md(job):
    with job:
        print(job)
    #print('Input file name given to runner is {}'.format(input_filename))
    #print('Output file name  is {}'.format(output_filename))


    #process=Popen("mpirun -n {} ~/test-cp2k/cp2k/exe/Linux-x86-64-intel/cp2k.popt -i {} -o {}".format(np,input_filename,output_filename),shell=True, universal_newlines=True,stdin=PIPE, stdout=PIPE, stderr=PIPE )
        process=Popen("lmp < control.inp", shell=True, universal_newlines=True,stdin=PIPE, stdout=PIPE, stderr=PIPE )
        output, error = process.communicate();
        print (output, error);
   # return output,error




if __name__ == "__main__":
    Project().main()

