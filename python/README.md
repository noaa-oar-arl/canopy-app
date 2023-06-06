## Getting started

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/noaa-oar-arl/canopy-app/HEAD?labpath=python%2Fexamples.ipynb)

Try it out on Binder :point_up: or follow the steps below :point_down: to run locally.

---

Check that your active Python environment has the dependencies listed in [`environment.yml`](./environment.yml).
Otherwise, you can create a new environment using:

    conda env create -f environment.yml

And then activate it using:

    conda activate canopy-app

The primary function is `run()`, which runs canopy-app and returns an xarray Dataset.
By default, it runs the [default case](../input/namelist.canopy).
```python
from canopy_app import run

ds = run()
```

You can override default namelist options by passing a dictionary to `run()`.
```python
# Point setup
ds = run(
    config={
        "filenames": {"file_vars": "../input/input_variables_point.txt"},
        "userdefs": {"infmt_opt": 1, "nlat": 1, "nlon": 1},
    },
)
```

There are also helper functions for running sets of experiments with different namelist options.
```python
# Multiple point setup cases
from canopy_app import config_cases, run_config_sens

cases = config_cases(
    file_vars="../input/input_variables_point.txt",
    infmt_opt=1,
    nlat=1,
    nlon=1,
    z0ghc=[0.001, 0.01],
    lambdars=[1.0, 1.25],
    product=True,  # 4 cases if true, 2 if false
)
ds = run_config_sens(cases)
```
:point_up: We still get a single output dataset, but it has a `case` dimension.

:point_down: You can also download and generate global gridded canopy-app inputs using python.

- Edit python script (`global_data_process.py`)

```
vim global_data_process.py
```

- Change user settings

```
'''User Options'''
path     = '/scratch/pcampbe8/canopy-app/input'    # work directory
year     = 2022   # year
month    = 7      # month
day      = 1      # day
houri    = 12     # gfs initialization hour in UTC (caution currently GFS files are initialized at 12 UTC only -- do not change)
hour     = 00     # gfs forecast hour in UTC
ref_lev  = 10     # reference height (m, a.g.l.)
frp_src  = 1      # frp data source (0: local source; 1: check local source first, switch to climatological file if no available data ; 2: 12 month climatology)
```

- Activate Python canopy-app environment:


```
    conda activate canopy-app
```

- Run python script

```
python3 global_data_process.py
```
