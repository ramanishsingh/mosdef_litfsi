import flow
from flow import FlowProject, directives
import warnings
import shutil
import os
from get_sol_il_xml import GetSolv, GetIL, Get_ff_path
import signac
import environment_for_rahman
from mtools.gromacs.gromacs import make_comtrj
import mdtraj as md
import numpy as np
from mtools.post_process import calc_msd

warnings.filterwarnings("ignore", category=DeprecationWarning)

has_md_files = "system.gro"
nvt_file = "nvt.gro"
sample_file = "sample.gro"
unwrapped_file = 'sample_unwrapped.xtc'

class Project(FlowProject):
    pass


@Project.label
def has_setter(job):
    return job.isfile("setter.py")


@Project.label
def has_structure(job):
    return job.isfile("init.xyz")



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
        shutil.copyfile(Project().root_directory() + "/init_structures/"+"struc_10m_small_{}K_{}.xyz".format(T,seed), job.workspace() + "/init.xyz")



@Project.label
def initialized(job):
    return job.isfile(has_md_files)


@Project.operation
@Project.post.isfile(has_md_files)
def initialize(job):

    import os
    import glob
    import numpy as np
    import unyt as u
    import mbuild as mb
    import ilforcefields.ilforcefields as ilff
    from foyer import Forcefield
    
    ### important setting
    case = '20m'
    
    
    cation = GetIL("Li")
    cation.name = "Li"
    anion = GetIL("TF2")
    anion.name = "TF2"
    sol = GetIL("wat")
    sol.name = "wat"
    
    
    if case == '10m':
        #### molecule number setting for 10m, 
        n_list = [9, 9, 50] #[Li, TF2, wat]
    if case == '20m':
        #### molecule number setting for 10m, 
        n_list = [14, 14, 39] #[Li, TF2, wat]
        
    print('n_list is ', n_list)
    
    packing_box = mb.Box([7,7,7]) ## note: only for eaiser to fill box, not final box size
    

        
    print("start to fill box")
    box_system = mb.fill_box(
                                compound = [cation, anion, sol],
                                n_compounds = [
                                            int(n_list[0]),
                                            int(n_list[1]),
                                            int(n_list[2])
                                            ],
                                box = packing_box,
                                edge=0.02
                                )
    cation_cmp = mb.Compound()
    anion_cmp = mb.Compound()
    sol_cmp = mb.Compound()
    
    for child in box_system.children:
        if child.name == 'Li':
            cation_cmp.add(mb.clone(child))
        elif child.name == 'TF2':
            anion_cmp.add(mb.clone(child))
        elif child.name == 'wat':
            sol_cmp.add(mb.clone(child)) 
    print('start to apply force field') 
    opls_li = Get_ff_path('opls_ions')
    opls_li = Forcefield(opls_li)
    cationPM = opls_li.apply(cation_cmp, residues = 'Li')
    clp = ilff.load_LOPES()
    anionPM = clp.apply(anion_cmp, residues = 'TF2')
    spce = Get_ff_path('spce')
    spce = Forcefield(spce)
    solPM = spce.apply(sol_cmp, residues = 'wat')
    structure = cationPM + anionPM + solPM
    
    print('start to scale partial charge to 0.8')
    for atom in structure.atoms:
        if atom.residue.name in ['Li','TF2']:
            atom.charge *= 0.8
            
    structure.combining_rule = 'geometric'
    
    init_pos_struc = mb.load('init_structures/struc_{}_small_{}K_{}.xyz'.format(case, int(job.statepoint()['T']), job.statepoint()['Seed']))
    init_pos_struc = init_pos_struc.to_parmed()
    pos = init_pos_struc.coordinates
    
    ### copy positions from cp2k to current md structure
    for i, atom in enumerate(structure.atoms):
        atom.xx = pos[i][0]
        atom.xy = pos[i][1]
        atom.xz = pos[i][2]
    
    if case == '10m':
        if job.statepoint()['T'] == 298:   
            box_size = [1.54361, 1.53940,1.54380]
            structure.box[0] = box_size[0] * 10
            structure.box[1] = box_size[1] * 10
            structure.box[2] = box_size[2] * 10
        
        if job.statepoint()['T'] == 373: 
            box_size = [1.57932,1.57932,1.57932]
            structure.box[0] = box_size[0] * 10
            structure.box[1] = box_size[1] * 10
            structure.box[2] = box_size[2] * 10
    
    if case == '20m':
        if job.statepoint()['T'] == 298:   
            box_size = [1.65843,1.65843,1.65843]
            structure.box[0] = box_size[0] * 10
            structure.box[1] = box_size[1] * 10
            structure.box[2] = box_size[2] * 10
        
        if job.statepoint()['T'] == 373: 
            box_size = [1.69073, 1.69207, 1.69298]
            structure.box[0] = box_size[0] * 10
            structure.box[1] = box_size[1] * 10
            structure.box[2] = box_size[2] * 10
    
    structure.save(os.path.join(job.workspace(), "system.gro"), combine ="all", overwrite=True)
    structure.save(os.path.join(job.workspace(), "system.top"), combine ="all", overwrite=True)
    


@Project.label
def nvt_done(job):
    return job.isfile(nvt_file)

@Project.operation
@Project.pre.isfile(has_md_files)
@Project.post.isfile(nvt_file)
@flow.cmd
def nvt(job):
    return _gromacs_str("nvt", "init", "init", job)

@Project.label
def sampled(job):
    return job.isfile(sample_file)

@Project.operation
@Project.pre.isfile(nvt_file)
@Project.post.isfile(sample_file)
@flow.cmd
def sample(job):
    return _gromacs_str("sample", "nvt", "init", job)


@Project.operation
@Project.pre.isfile(has_md_files)
@Project.post.isfile(sample_file)
@flow.cmd
def sample_total(job):
    return _gromacs_str_total('nvt','sample', job.statepoint()["T"], job)


@Project.label
def prepared(job):
    return job.isfile(unwrapped_file)

@Project.operation
@Project.pre.isfile(sample_file)
@Project.post.isfile(unwrapped_file)
def prepare(job):
    xtc_file = os.path.join(job.workspace(), 'sample.xtc')
    gro_file = os.path.join(job.workspace(), 'sample.gro')
    tpr_file = os.path.join(job.workspace(), 'sample.tpr')
    if os.path.isfile(xtc_file) and os.path.isfile(gro_file):
        unwrapped_trj = os.path.join(job.workspace(),
        'sample_unwrapped.xtc')
        os.system('echo 0 | gmx trjconv -f {0} -o {1} -s {2} -pbc nojump'.format(xtc_file, unwrapped_trj, tpr_file))
        res_trj = os.path.join(job.ws, 'sample_res.xtc')
        com_trj = os.path.join(job.ws, 'sample_com.xtc')
        unwrapped_com_trj = os.path.join(job.ws,'sample_com_unwrapped.xtc')
        os.system('echo 0 | gmx trjconv -f {0} -o {1} -s {2} -pbc res'.format(
                xtc_file, res_trj, tpr_file))
        trj = md.load(res_trj, top=gro_file)
        comtrj = make_comtrj(trj)
        comtrj.save_xtc(com_trj)
        comtrj[-1].save_gro(os.path.join(job.workspace(),
             'com.gro'))
        print('made comtrj ...')
        os.system('gmx trjconv -f {0} -o {1} -pbc nojump'.format(
                com_trj, unwrapped_com_trj))

@Project.label
def msd_done(job):
    return job.isfile(msd_file)

msd_file = 'diffusivity_and_msd_done.txt_1'
@Project.operation
@Project.pre.isfile(unwrapped_file)
@Project.post.isfile(msd_file)
def run_msd(job):
    import matplotlib.pyplot as plt
    print('Loading trj {}'.format(job))
    top_file = os.path.join(job.workspace(), 'com.gro')
    trj_file = os.path.join(job.workspace(),
            'sample_com_unwrapped.xtc')
    trj = md.load(trj_file, top=top_file)


    selections = {'all' : trj.top.select('all'),
                'ion' : trj.top.select('resname li tfsi'),
                'cation': trj.top.select("resname li"),
                'anion': trj.top.select("resname tfsi"),
                'solvent': trj.top.select("resname wat"),
                }

    for mol, indices in selections.items():
        print('\tConsidering {}'.format(mol))
        if indices.size == 0:
            print('{} does not exist in this statepoint'.format(mol))
            continue
        print(mol)
        sliced = trj.atom_slice(indices)
        D, MSD = _run_overall(sliced, mol)
        job.document['D_' + mol + '_overall_original'] = D
        
        np.savetxt('msd_data/msd_{}K_{}.txt'.format(job.statepoint()["T"], job.statepoint()["Seed"]),np.transpose(np.vstack([trj.time, MSD])),header='# Time (ps)\tMSD (nm^2)')
        fig, ax = plt.subplots()
        ax.plot(trj.time, MSD)
        ax.set_xlabel('Simulation time (ps)')
        ax.set_ylabel('MSD (nm^2)')
        fig.savefig('msd_data/msd_{}K_{}.jpg'.format(job.statepoint()["T"], job.statepoint()["Seed"]))

    np.savetxt(os.path.join(job.workspace(), 'diffusivity_and_msd_done.txt_1'), [1,1]) 
    
       
def _run_overall(trj, mol):
    D, MSD, x_fit, y_fit = calc_msd(trj)
    return D, MSD

def workspace_command(cmd):
    """Simple command to always go to the workspace directory"""
    return " && ".join(
        [
            "cd {job.ws}",
            cmd if not isinstance(cmd, list) else " && ".join(cmd),
            "cd ..",
        ]
    )
    
def _gromacs_str(op_name, gro_name, sys_name, job):
    """Helper function, returns grompp command string for operation """
    if op_name == 'em':
        mdp = signac.get_project().fn('util/mdp_files/{}.mdp'.format(op_name))
        cmd = ('gmx grompp -f {mdp} -c system.gro -p system.top -o {op}.tpr --maxwarn 1 && gmx mdrun -deffnm {op} -ntmpi 1')
    else:
        mdp = signac.get_project().fn('util/mdp_files/{}.mdp'.format(op_name))
        cmd = ('gmx grompp -f {mdp} -c {gro}.gro -p system.top -o {op}.tpr --maxwarn 1 && gmx mdrun -deffnm {op} -ntmpi 1')
    return workspace_command(cmd.format(mdp=mdp,op=op_name, gro=gro_name, sys=sys_name))


def _gromacs_str_total(op_name_1, op_name_2, T, job):
    """Helper function, returns grompp command string for operation """
    mdp1 = signac.get_project().fn('util/mdp_files/{}-{}.mdp'.format(op_name_1, int(T)))
    mdp2 = signac.get_project().fn('util/mdp_files/{}-{}.mdp'.format(op_name_2, int(T)))
    cmd = (
            'gmx grompp -f {mdp1} -c {gro1}.gro -p system.top -o {op1}.tpr --maxwarn 1 '\
            '&& gmx mdrun -deffnm {op1} -ntmpi 1 '\
            '&& gmx grompp -f {mdp2} -c {gro2}.gro -p system.top -o {op2}.tpr --maxwarn 1 '\
            '&& gmx mdrun -deffnm {op2} -ntmpi 1'
            )
    return workspace_command(cmd.format(mdp1=mdp1, op1 = op_name_1, gro1 = 'system',
                                        mdp2=mdp2, op2 = op_name_2, gro2 = op_name_1,
                                        ))

if __name__ == "__main__":
    Project().main()

