import glob, os, sys
import numpy as np

import xarray as xr
import pandas as pd
import pyarrow.parquet as pq

import scipy.ndimage.measurements as measure
import scipy.ndimage.morphology as morph

from pathlib import Path
from tqdm import tqdm

# Define conditional cloud field
cp = 1004.0  # Heat capacity at constant pressure for dry air [J kg^-1 K^-1]
Rv = 461.0  # Gas constant of water vapor [J kg^-1 K^-1]
Rd = 287.0  # Gas constant of dry air [J kg^-1 K^-1]
p_0 = 100000.0
lam = Rv / Rd - 1.0


def theta_v(p, T, qv, qn, qp):
    return T * (1.0 + lam * qv - qn - qp) * (p_0 / p) ** (Rd / cp)


def sample_cond(ds):
    """
    Define conditional fields

    Return
    ------
    c0_fld: cloud field (QN > 0)
    c1_fld: core field (QN > 0, W > 0, B > 0)
    c2_fld: plume (tracer-dependant)
    """
    th_v = theta_v(
        ds["p"][:] * 100,
        ds["TABS"][:],
        ds["QV"][:] / 1e3,
        ds["QN"][:] / 1e3,
        ds["QP"][:] / 1e3,
    )

    buoy = th_v > np.mean(th_v, axis=(1, 2))

    c0_fld = ds["QN"] > 0
    c1_fld = buoy & c0_fld

    # Define plume based on tracer fields (computationally intensive)
    # tr_field = ds['TR01'][:]
    # tr_ave = np.nanmean(tr_field, axis=(1, 2))
    # tr_std = np.std(tr_field, axis=(1, 2))
    # tr_min = .05 * np.cumsum(tr_std) / (np.arange(len(tr_std))+1)

    # c2_fld = (tr_field > \
    #             np.max(np.array([tr_ave + tr_std, tr_min]), 0)[:, None, None])

    return c0_fld, c1_fld  # , c2_fld


def cluster_clouds(ds):
    """
    Clustering the cloud volumes/projections from Zarr dataset.

    Input
    -----
    ds: Zarr file

    Return
    ------
    Parquet files containing the coordinates
        coord: point coordinates (see lib/index.py)  
        cid: cloud id
        type: 
            - 0 for condensed region 
            - 1 for core region
            - 2 for plumes 

    """
    bin_st = morph.generate_binary_structure(3, 2)
    bin_st[0] = 0
    bin_st[-1] = 0

    c0_fld, c1_fld = sample_cond(ds)

    df = pd.DataFrame(columns=["coord", "cid", "type"])
    for item in [0, 1]:
        c_field = locals()[f"c{item}_fld"]

        c_label, n_features = measure.label(c_field, structure=bin_st)
        c_label = c_label.ravel()  # Sparse array

        # Extract indices
        c_index = np.arange(len(c_label))
        c_index = c_index[c_label > 0]
        c_label = c_label[c_label > 0]

        if item == 0:
            c_type = np.ones(len(c_label), dtype=int)
        elif item == 1:
            c_type = np.zeros(len(c_label), dtype=int)
        else:
            raise TypeError

        df_ = pd.DataFrame.from_dict(
            {"coord": c_index, "cid": c_label, "type": c_type,}
        )
        df = pd.concat([df, df_])
    return df


if __name__ == "__main__":
    # TODO: Use ModelConfig class for model parameters
    case_name = "BOMEX_SWAMP"

    src = Path(f"/Howard16TB/data/loh/{case_name}/variables")
    dest = "/users/loh/nodeSSD/temp"

    fl = sorted(src.glob("*"))
    for time, item in enumerate(tqdm(fl)):
        fname = item.as_posix()
        if item.suffix == ".zarr":
            with xr.open_zarr(fname) as ds:
                try:
                    ds = ds.squeeze("time")
                except:
                    pass

            df = cluster_clouds(ds)
            print(df)
            raise
        else:
            print("Error: File type not recognized.")
            raise Exception
