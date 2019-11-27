import xarray
from netCDF4 import Dataset
from pathlib import Path
import context
from datetime import datetime
from pytz import utc
import matplotlib
from matplotlib import pyplot as plt
import pprint
import xarray as xr
import pdb
matplotlib.use("Agg")




if __name__ == "__main__":
    zarr_in = list(context.data_dir.glob("*zarr"))[0]
    ds_zarr = xr.open_zarr(str(zarr_in))
    pdb.set_trace()
