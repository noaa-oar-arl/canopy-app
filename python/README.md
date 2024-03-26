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
        "filenames": {"file_vars": "../input/point_file_20220701.sfcf000.txt"},
        "userdefs": {"infmt_opt": 1, "nlat": 1, "nlon": 1},
    },
)
```

There are also helper functions for running sets of experiments with different namelist options.
```python
# Multiple point setup cases
from canopy_app import config_cases, run_config_sens

cases = config_cases(
    file_vars="../input/point_file_20220701.sfcf000.txt",
    infmt_opt=1,
    ntime=1,
    nlat=1,
    nlon=1,
    z0ghc=[0.001, 0.01],
    lambdars=[1.0, 1.25],
    product=True,  # 4 cases if true, 2 if false
)
ds = run_config_sens(cases)
```
:point_up: We still get a single output dataset, but it has a `case` dimension.
