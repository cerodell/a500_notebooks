# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
''''
demo
'''
import xarray
from pathlib import Path
import xarray as xr
import re
import os

timestep = re.compile(r".*_(\d+).zarr")


def sort_step(filename):
    out = str(filename)
    the_match = timestep.match(out)
    if the_match:
        the_step = int(the_match.group(1))
    else:
        the_step = None
    return the_step


if __name__ == "__main__":
    data_dir = Path("/scratch/paustin/paustin/a500_data/small_bomex")
    zarr_in = list(data_dir.glob("*zarr"))
    zarr_in.sort(key=sort_step)
    print(zarr_in[:10])
    zarr_list = []
    the_time = 0
    #
    # open each zarr file as an xarray dataset
    # correcting the timestamp
    #
    for ds_zarr in zarr_in:
        ds = xr.open_zarr(str(ds_zarr))
        #
        # wrote incorrect times in original files, fix here
        #
        ds['time'] = the_time
        zarr_list.append(ds)
        the_time += 60
    #
    # make a virtual dataset with time as the outer dimension
    #
    zarr_time_ds = xr.combine_nested(zarr_list, 'time')
    #
    # compute the mean and perturbation for timestep 0 and
    # write out as a new zarr file
    #
    time_step = 0
    print(f"finding perturbation for  {zarr_in[time_step]}")
    temp = zarr_time_ds['TABS']
    wvel = zarr_time_ds['W']
    tr01 = zarr_time_ds['TR01']
    mean_temp = temp[time_step, :, :, :].mean(dim=('x', 'y'))
    mean_w = wvel[time_step, :, :, :].mean(dim=('x', 'y'))
    mean_tr = tr01[time_step, :, :, :].mean(dim=('x', 'y'))
    w_prime = wvel - mean_w
    T_prime = temp - mean_temp
    tr_prime = tr01 - mean_tr
    #
    # put these variables in a dictionary so they can be placed
    # in a new dataset
    #
    varnames = [
        'mean_temp', 'mean_w', 'mean_tr', 'w_prime', 'T_prime', 'tr_prime'
    ]
    ds_vec = [mean_temp, mean_w, mean_tr, w_prime, T_prime, tr_prime]
    new_ds = xr.Dataset(dict(zip(varnames, ds_vec)))
    #
    # keep this small for illustration purposes
    #
    new_ds_slim = new_ds.isel(x=slice(0, 40), y=slice(0, 50))
    #
    # make the timestep from the original file part of the name.
    #
    parts = str(zarr_in[time_step]).split("_")
    user_name = os.environ['USER']
    out_dir = Path(f"/scratch/paustin/{user_name}/flux_bomex")
    out_dir.mkdir(parents=True, exist_ok=True)
    full_path = out_dir / f"flux_{parts[-1]}"
    #
    # the next step is important!!!  if you write out the ds
    # without first computing, it will write the entire multi-timestep file
    #
    new_ds_slim.compute()
    new_ds_slim.to_zarr(full_path, "w")
    print(f"wrote {full_path}")
