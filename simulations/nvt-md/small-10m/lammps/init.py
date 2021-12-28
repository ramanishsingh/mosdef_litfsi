import signac
import numpy as np
import unyt as u

def init_project():

    # Initialize project
    project = signac.init_project("small-10m")

    # Define temperature
    temperatures = [298, 373]* u.K
    lengths=np.array([[1.54361,1.53940,1.54380], [1.57932,1.57932,1.57932]])*u.nm
    seeds=[1,2,3,4]
    i=0
    for temperature in temperatures:
        for seed in seeds:
            # Define the state point
            state_point = {
                "T": float(temperature.in_units(u.K).value),
                "L": lengths[i].in_units(u.Angstrom).value,
                "Seed" : seed
            }

            job = project.open_job(state_point)
            job.init()
        i+=1


if __name__ == "__main__":
    init_project()

