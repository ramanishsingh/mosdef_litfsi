import numpy as np

processes = ['nvt', 'npt', 'equil', 'sample']

temperatures = [260, 280, 300]

original = 'ref_t               = 300\n'

for proc in processes:
    for temp in temperatures:
        new = 'ref_t               = {}'.format(temp)
        with open('{}.mdp'.format(proc), 'r') as f_in, open('{}-{}.mdp'.format(proc,
            temp), 'w') as f_out:
            for line_in in f_in:
                if line_in == original:
                    f_out.write(new)
                else:
                    f_out.write(line_in)