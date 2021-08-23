import signac
import mbuild
import numpy as np
import pandas as pd
import unyt as u
import mdtraj as md

from mosdef_slitpore.utils.cassandra_helpers import spce_water
from mosdef_slitpore.analysis import compute_density
from mosdef_slitpore.analysis import compute_s


def main():

    project = signac.get_project("../")


    for nwater, group in project.groupby("nwater"):
        # New dataframe to save the results
        df = pd.DataFrame()

        # Create pore system
        pore_width = 2.0 * u.nm

        filled_pore = mbuild.recipes.GraphenePoreSolvent(
            pore_length=3.0,
            pore_depth=3.0,
            pore_width=pore_width.to_value("nm"),
            n_sheets=3,
            slit_pore_dim=2,
            x_bulk=0,
            solvent=spce_water,
            n_solvent=nwater,
        )
        # Translate to centered at 0,0,0 and make box larger in z
        box_center = filled_pore.periodicity/2.0
        filled_pore.translate(-box_center)
        filled_pore.periodicity[2] = 6.0

        xy_area = filled_pore.periodicity[0] * filled_pore.periodicity[1]
        top = filled_pore.to_trajectory(residues=["RES", "SOL"])

        # Load all trajectories and combine
        for job in group:
            run = job.sp.run
            # Load in full trajectory
            full_traj = md.load(job.fn("prod.nvt.out.xyz"), top=top)
            # Add unit cell information
            full_traj = md.Trajectory(
                full_traj.xyz,
                full_traj.top,
                unitcell_lengths = np.tile(filled_pore.periodicity, (full_traj.n_frames,1)),
                unitcell_angles = np.tile([90.,90.,90.], (full_traj.n_frames,1)),
            )
            # Keep only water
            slice_water = full_traj.top.select("water and name O H")
            traj = full_traj.atom_slice(slice_water)
            slice_ow = traj.top.select("name O")
            slice_hw = traj.top.select("name H")
            traj_ow = traj.atom_slice(slice_ow)
            traj_hw = traj.atom_slice(slice_hw)

            # Compute the density
            bin_centers_ow, density_ow = compute_density(traj_ow, xy_area, bin_width=0.005)
            bin_centers_hw, density_hw = compute_density(traj_hw, xy_area, bin_width=0.005)

            # Compute the s order parameter
            bin_centers_s, s_results = compute_s(traj, bin_width=0.005)

            assert np.allclose(bin_centers_ow, bin_centers_hw)
            assert np.allclose(bin_centers_ow, bin_centers_s)

            # Save results
            tmp_df = pd.DataFrame()
            tmp_df[f"run"] = run
            tmp_df[f"z-loc_nm"] = bin_centers_ow
            tmp_df[f"density-ow_nm^-3"] = density_ow
            tmp_df[f"density-hw_nm^-3"] = density_hw
            tmp_df[f"s_value"] = s_results
            df = df.append(tmp_df)

        # Compute mean/stdev and save to file
        means = df.groupby("z-loc_nm").mean().drop(columns=["run"]).add_suffix("_mean")
        stds = df.groupby("z-loc_nm").std().drop(columns=["run"]).add_suffix("_std")
        combined = means.merge(stds, on="z-loc_nm")
        combined.to_csv(f"results_{nwater}-water.csv")


if __name__ == "__main__":
    main()
