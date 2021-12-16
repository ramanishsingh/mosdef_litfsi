import flow
from flow import FlowProject, directives
import warnings
import shutil
import os

warnings.filterwarnings("ignore", category=DeprecationWarning)


class Project(FlowProject):
    pass


@Project.label
def has_setter(job):
    return job.isfile("setter.py")


@Project.label
def has_structure(job):
    return job.isfile("init.xyz")


@Project.label
def has_md_files(job):
    "Verify that the files required for MD simulation are there"

    try:
        return job.isfile(job.doc.input_filename)
    except:
        return False


@Project.label
def has_restart_file(job):
    "Verify that the restart file required for restarting the MD simulation is there"

    try:
        return job.isfile(job.doc.restart_filename)
    except:
        return False



@Project.label
def md_completed(job):
    "Verify that the md simulation has completed"

    try:
        return job.isfile(job.doc.output_filename)
    except:
        return False


@Project.operation
@Project.post(has_setter)
def copy_setter(job):
    with job:
        print(os.listdir())
        shutil.copyfile(
            Project().root_directory() + "/setter.py", job.workspace() + "/setter.py"
        )

@Project.operation
@Project.post(has_structure)
def copy_structure(job):
    with job:
        seed=int(job.sp.Seed)
        T=int(job.sp.T)
        shutil.copyfile(Project().root_directory() + "/init_structures/"+"struc_20m_small_{}K_{}.xyz".format(T,seed), job.workspace() + "/init.xyz")


@Project.operation
@Project.post(has_md_files)
def md_files(job):
    import mbuild as mb

    with job:
        import os
        import glob
        import numpy as np
        import unyt as u
        import mbuild as mb
        from cp2kmdpy.molecule_optimization import (
            Molecule_optimization,
        )  # for single molecule optimization
        from cp2kmdpy.md import MD  # for running MD
        from cp2kmdpy import runners
        import setter

        temperature = job.sp.T * u.K
        functional= job.sp.Functional
        length=job.sp.L
        print(length)
        seed=job.sp.Seed
        molecule = mb.load("init.xyz")
        box = mb.box.Box(lengths=[length[0]/10, length[1]/10, length[2]/10])
        q = MD(
            molecules=[molecule],
            box=box,
            cutoff=600,
            functional=functional,
            basis_set={
                "C": "DZVP-MOLOPT-SR-GTH",
                "H": "DZVP-MOLOPT-SR-GTH",
                "O": "DZVP-MOLOPT-SR-GTH",
                "Li": "DZVP-MOLOPT-SR-GTH",
                "F": "DZVP-MOLOPT-SR-GTH",
                "S": "DZVP-MOLOPT-SR-GTH",
                "N": "DZVP-MOLOPT-SR-GTH",
            },
            periodicity="XYZ",
            n_molecules=[1],
            traj_type="XYZ",
            seed=seed,
            project_name="litfsi",
            initial_coordinate_filename="init.xyz"
        )
        q.temperature = temperature
        q.ensemble = "NVT"
        q.simulation_time = 200 * u.ps
        # Initializing q
        q.md_initialization()

        # generating input files
        setter.md_files(q)
        job.doc.input_filename = q.input_filename
        job.doc.output_filename = q.output_filename
        job.doc.restart_filename = q.project_name + "-1.restart"


@Project.operation
@Project.pre(has_md_files)
@Project.post(has_restart_file)
@flow.directives(np=72)
def run_md(job):
    from cp2kmdpy import runners_mpi
    import os

    with job:

        a = runners_mpi.run_md(job.doc.input_filename, job.doc.output_filename, 72)
        print(a)


@Project.operation
@Project.pre(has_restart_file)
@flow.directives(np=72)
def restart_md(job):
    from cp2kmdpy import runners_mpi
    import os

    with job:

        a = runners_mpi.run_md(job.doc.restart_filename, job.doc.output_filename, 72)
        print(a)


if __name__ == "__main__":
    Project().main()
