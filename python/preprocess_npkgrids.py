"""
Preprocess NPKGRIDS fertilizer data for canopy-app soil NO emissions (BDSNP).

Downloads (or reads local) NPKGRIDS v1.08 NetCDF files (173 crop-specific),
sums nitrogen application rates (Nrate, kg-N/ha) across all crops per grid
cell, and regrids from native 0.05 deg to the GFS ~13 km grid using
monet.remap_nearest().

Usage:
    python preprocess_npkgrids.py <gfs_reference_file> [npkgrids_dir]

Arguments:
    gfs_reference_file : Path to any GFS meteorological NetCDF file (defines
                         the target 1536x3072 grid).
    npkgrids_dir       : Directory containing extracted NPKGRIDSv1.08_*.nc
                         files.  Default: ./input/npkgrids

Output:
    ./input/nitrogen_input.nc   (total N application rate on GFS grid)

Reference:
    Nguyen et al. (2024), NPKGRIDS: a global georeferenced dataset of N,
    P2O5, and K2O fertilizer application rates for 173 crops.
    https://doi.org/10.1038/s41597-024-04030-4
    Data: https://doi.org/10.6084/m9.figshare.24616050

Author: Quazi Rasool (CIRES/NOAA CSL)
"""

import glob
import os
import sys
from datetime import datetime

import monet  # noqa: F401
import numpy as np
import xarray as xr
from netCDF4 import Dataset

# ----------------------------- User arguments ------------------------------ #
if len(sys.argv) < 2:
    print("Usage: python preprocess_npkgrids.py <gfs_reference_file> [npkgrids_dir]")
    print("  gfs_reference_file : any GFS .nc file (for target grid)")
    print("  npkgrids_dir       : folder with NPKGRIDSv1.08_*.nc (default: ./input/npkgrids)")
    sys.exit(1)

f_gfs = sys.argv[1]
npk_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(".", "input", "npkgrids")
f_output = os.path.join(".", "input", "nitrogen_input.nc")

fill_value = 9.99e20

starttime = datetime.now()
print("------------------------------------")
print("---- NPKGRIDS Preprocessing")
print("---- Start time:", starttime.strftime("%Y/%m/%d %H:%M:%S"))
print("------------------------------------")

# ----------------------------- Locate files -------------------------------- #
nc_files = sorted(glob.glob(os.path.join(npk_dir, "NPKGRIDSv1.08_*.nc")))
if not nc_files:
    print(f"ERROR: No NPKGRIDSv1.08_*.nc files found in {npk_dir}")
    print("Download from https://doi.org/10.6084/m9.figshare.24616050")
    print("and extract NPKGRIDSv1.08_NC.zip into the directory above.")
    sys.exit(1)

print(f"---- Found {len(nc_files)} crop files in {npk_dir}")

# -------------------- Sum Nrate across all crops --------------------------- #
print("---- Summing Nrate across all crops at native 0.05 deg ...")

total_n = None
n_crops_counted = 0

for fpath in nc_files:
    crop_name = os.path.basename(fpath).replace("NPKGRIDSv1.08_", "").replace(".nc", "")
    ds = xr.open_dataset(fpath)

    if "Nrate" not in ds.data_vars:
        print(f"  SKIP {crop_name}: no Nrate variable")
        ds.close()
        continue

    nrate = ds["Nrate"].values.copy()
    ds.close()

    # Squeeze extra dimensions (some files may have time=1)
    nrate = np.squeeze(nrate)

    # Ocean / water cells are marked as -1; treat as zero
    nrate[nrate < 0] = 0.0
    # Also mask NaN
    nrate[np.isnan(nrate)] = 0.0

    if total_n is None:
        total_n = nrate.astype(np.float64)
    else:
        total_n += nrate.astype(np.float64)

    n_crops_counted += 1

print(f"---- Summed {n_crops_counted} crops.  Max total N = {total_n.max():.2f} kg-N/ha")

# -------------------- Build source xarray for monet ----------------------- #
# Re-open one crop file to get native coordinate arrays
ds0 = xr.open_dataset(nc_files[0])

# Detect coordinate names (lat/lon or latitude/longitude)
coord_names = list(ds0.coords) + list(ds0.dims)
lat_name = next((c for c in coord_names if c.lower() in ("lat", "latitude")), None)
lon_name = next((c for c in coord_names if c.lower() in ("lon", "longitude")), None)

if lat_name is None or lon_name is None:
    # Fallback: construct from 0.05 deg grid
    print("---- WARNING: Could not detect lat/lon coords; constructing 0.05-deg grid")
    nlat_src, nlon_src = total_n.shape
    lat_src = np.linspace(-90 + 0.025, 90 - 0.025, nlat_src)
    lon_src = np.linspace(-180 + 0.025, 180 - 0.025, nlon_src)
else:
    lat_src = ds0[lat_name].values
    lon_src = ds0[lon_name].values

ds0.close()

# Build a 2-D lat/lon mesh and wrap into an xarray DataArray
lon2d, lat2d = np.meshgrid(lon_src, lat_src)
src_da = xr.DataArray(
    total_n,
    dims=["y", "x"],
    coords={
        "latitude": (["y", "x"], lat2d),
        "longitude": (["y", "x"], lon2d),
    },
)

print(f"---- Source grid shape: {total_n.shape}  ({len(lat_src)} lat x {len(lon_src)} lon)")

# -------------------- Open GFS reference and regrid ------------------------ #
print(f"---- Reading GFS reference grid from {f_gfs} ...")
basefile = xr.open_dataset(f_gfs)
basefile = basefile.set_coords(["lat", "lon"]).rename(
    {"grid_xt": "x", "grid_yt": "y", "lat": "latitude", "lon": "longitude"}
)

nlat_gfs = len(basefile["y"].data)
nlon_gfs = len(basefile["x"].data)
print(f"---- Target GFS grid: {nlat_gfs} x {nlon_gfs}")

# Find a 2-D reference variable for monet remap target grid
# (global_data_process.py uses basefile["zc"]; fall back to any 2-D var)
ref_var = None
for vn in ["zc", "tmpsfc", "tmp2m", "pressfc"]:
    if vn in basefile.data_vars:
        ref_var = vn
        break
if ref_var is None:
    # pick first 2-D (or 3-D with time) data variable
    for vn in basefile.data_vars:
        ndim = len(basefile[vn].dims)
        if ndim in (2, 3):
            ref_var = vn
            break
if ref_var is None:
    print("ERROR: Could not find a suitable 2-D reference variable in GFS file")
    sys.exit(1)

ref_da = basefile[ref_var]
if len(ref_da.dims) == 3:  # (time, y, x) → take first time slice
    ref_da = ref_da[0, :, :]

print(f"---- Using '{ref_var}' as reference grid variable")
print("---- Regridding with monet.remap_nearest() ...")
ninput_gfs = ref_da.monet.remap_nearest(src_da).data

# Clean up regridded result
ninput_gfs[np.isnan(ninput_gfs)] = 0.0
ninput_gfs[ninput_gfs < 0] = 0.0

print(f"---- Regridded N-input: min={ninput_gfs.min():.2f}, max={ninput_gfs.max():.2f}, "
      f"mean={ninput_gfs.mean():.2f} kg-N/ha")

basefile.close()

# -------------------- Write output NetCDF ---------------------------------- #
print(f"---- Writing output to {f_output} ...")

# Re-open GFS file to copy coordinate variables
gfs_ds = Dataset(f_gfs, "r")
out_ds = Dataset(f_output, "w", format="NETCDF4")

# Copy dimensions
out_ds.createDimension("grid_yt", nlat_gfs)
out_ds.createDimension("grid_xt", nlon_gfs)

# Copy lat/lon coordinate variables from GFS
for vname in ["grid_xt", "grid_yt", "lat", "lon"]:
    src_var = gfs_ds.variables[vname]
    out_var = out_ds.createVariable(vname, src_var.dtype, src_var.dimensions)
    out_var[:] = src_var[:]
    for attr in src_var.ncattrs():
        out_var.setncattr(attr, src_var.getncattr(attr))

# Write total N-input variable
var = out_ds.createVariable(
    "ninput", "f4", ("grid_yt", "grid_xt"), fill_value=fill_value
)
var.long_name = "Total nitrogen fertilizer application rate (all crops)"
var.units = "kg-N ha-1"
var.missing_value = fill_value
var.source = "NPKGRIDS v1.08 (Nguyen et al., 2024)"
var.reference = "https://doi.org/10.1038/s41597-024-04030-4"
var[:] = ninput_gfs

# Global attributes
out_ds.title = "Nitrogen input for BDSNP soil NO emissions (canopy-app)"
out_ds.source = f"Preprocessed from NPKGRIDS v1.08 ({n_crops_counted} crops summed)"
out_ds.history = f"Created {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by preprocess_npkgrids.py"
out_ds.conventions = "CF-1.6"
out_ds.reference = "Nguyen et al. (2024), Sci Data 11:1179, doi:10.1038/s41597-024-04030-4"

gfs_ds.close()
out_ds.close()

endtime = datetime.now()
print("------------------------------------")
print(f"---- Output written: {f_output}")
print(f"---- Grid: {nlat_gfs} x {nlon_gfs} (GFS ~13 km)")
print(f"---- Process time: {endtime - starttime}")
print("---- NPKGRIDS preprocessing complete!")
print("------------------------------------")

