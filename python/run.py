"""
Run canopy-app using the namelist and load nc data.
"""
from __future__ import annotations

import contextlib
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Callable

import f90nml
import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).parent.absolute()
assert HERE.name == "python"
assert HERE.parent.name == "canopy-app"
REPO = HERE.parent

with open(REPO / "input" / "namelist.canopy") as f:
    DEFAULTS = f90nml.read(f)

del f

# Make input paths absolute in default config
for k, v in DEFAULTS["filenames"].items():
    p0 = Path(v)
    if not p0.is_absolute():
        p = REPO / p0
    else:
        p = p0
    DEFAULTS["filenames"][k] = p.as_posix()

del (k, v, p0, p)


@contextlib.contextmanager
def out_and_back(p: Path, *, finally_: Callable | None = None):
    """Context manager: change working directory to path `p` but return to original cwd after."""
    cwd = os.getcwd()
    os.chdir(p)
    try:
        yield
    finally:
        if finally_ is not None:
            finally_()
        os.chdir(cwd)


POINT_DEFAULT = pd.read_csv(REPO / "input" / "input_variables_point.txt", index_col=False)

_TXT_STEM_SUFFS = {
    "wind": "_output_canopy_wind",
    "waf": "_output_waf",
    "eddy": "_output_eddy_Kz",
    "phot": "_output_phot",
    "bio": "_output_bio",
}


def run(
    *,
    config: dict[str, Any] | None = None,
    case_dir: Path | None = None,
    cleanup: bool = True,
    verbose: bool = False,
) -> xr.Dataset:
    """Run canopy-app, returning xarray dataset.

    Parameters
    ----------
    config
        Namelist options as a dict.
        Defaults used if not provided.
        (Don't need to change input path, we set to absolute path automatically.)
    case_dir
        Directory for the case.
        A temporary directory is used if not provided.
        Replicates a bit of the structure (``input``, ``output`` subdirectories).
    cleanup
        Clean up (remove) the case directory after (successful) run.
    """
    # Ready the output directory
    if case_dir is None:
        import tempfile

        temp_dir = tempfile.TemporaryDirectory(prefix="canopy-app_")
        case_dir = Path(temp_dir.name)
    else:
        case_dir.mkdir(parents=True, exist_ok=True)

    # Create config
    full_config: f90nml.Namelist = DEFAULTS.copy()
    user_config = config or {}
    for section_name, sub_config in user_config.items():
        if section_name not in full_config:
            raise ValueError(f"Unexpected namelist section: {section_name!r}")
        for k, v in sub_config.items():
            assert not isinstance(v, dict)
            if k not in full_config[section_name]:
                raise ValueError(f"Unexpected namelist option: '{section_name}.{k}'")
            if verbose:
                print(f"'{section_name}.{k}' {full_config[section_name][k]} -> {v!r}")
            full_config[section_name][k] = v
    output_dir = case_dir / "output"
    output_dir.mkdir(exist_ok=True)
    ofp_stem = output_dir / "out"
    full_config["filenames"]["file_out"] = ofp_stem.relative_to(case_dir).as_posix()

    # Check input file
    ifp = Path(full_config["filenames"]["file_vars"])
    if not ifp.is_file():
        raise ValueError(f"Input file {ifp.as_posix()!r} does not exist")
    if not ifp.is_absolute():
        full_config["filenames"]["file_vars"] = ifp.resolve().as_posix()
    if ifp.suffix in {".nc", ".nc4", ".ncf"}:  # consistent with sr `canopy_check_input`
        nc_out = True
        assert full_config["userdefs"]["infmt_opt"] == 0
    elif ifp.suffix in {".txt"}:
        nc_out = False
        assert full_config["userdefs"]["infmt_opt"] == 1
    else:
        raise ValueError(f"Unexpected input file type: {ifp.suffix}")

    # Write namelist
    if verbose:
        print("Full namelist config:")
        print(full_config)
    config_dir = case_dir / "input"
    config_dir.mkdir(exist_ok=True)
    full_config.write(config_dir / "namelist.canopy", force=True)

    # Run canopy-app
    exe = REPO / "canopy"
    if not exe.is_file():
        raise RuntimeError("compile canopy-app first")
    cmd = [exe.as_posix()]
    with out_and_back(case_dir):
        subprocess.run(cmd, check=True, capture_output=not verbose)

    # Load nc
    if nc_out:
        ds0 = xr.open_dataset(ofp_stem.with_suffix(".nc"))
        ds = (
            ds0.rename_dims(grid_xt="x", grid_yt="y")
            .swap_dims(level="z")
            .assign(time=ds0.time)
            .set_coords(["time", "z", "lat", "lon"])
            .squeeze()
        )
    else:
        dfs: list[pd.DataFrame] = []
        ifcans = [
            k
            for k, v in full_config["userdefs"].items()
            if k.startswith("ifcan") and v is True
        ]
        if not ifcans:
            raise RuntimeError("None of the 'ifcan' switches are on, so no text output.")

        # Read txt outputs
        for ifcan in ifcans:
            which = ifcan[len("ifcan") :]
            stem_suff = _TXT_STEM_SUFFS.get(which)
            if stem_suff is None:
                print(
                    f"warning: skipping {which!r} ({ifcan}) output since stem suffix unknown."
                )
                continue
            df = read_txt(
                ofp_stem.with_name(f"{ofp_stem.name}{stem_suff}").with_suffix(".txt")
            )
            df.attrs.update(which=which)
            dfs.append(df)

        # Merge
        units: dict[str, str] = {}
        dss = []
        for df in dfs:
            if {"lat", "lon", "height"}.issubset(df.columns):
                ds_ = df.set_index(["height", "lat", "lon"]).to_xarray().squeeze()
            elif {"lat", "lon"}.issubset(df.columns):
                ds_ = df.set_index(["lat", "lon"]).to_xarray().squeeze()
            else:
                raise ValueError("Expected df to have columns 'lat', 'lon' [,'height'].")
            units.update(df.attrs["units"])
            for vn in ds_.data_vars:
                ds_[vn].attrs["units"] = units[vn]
                ds_[vn].attrs["group"] = df.attrs["which"]
            ds_.attrs.update(
                {k: v for k, v in df.attrs.items() if k not in {"which", "units"}}
            )
            dss.append(ds_)
        ds = xr.merge(dss, combine_attrs="no_conflicts")
        ds = ds.rename(height="z")
        if {"lat", "lon"}.issubset(ds.dims):
            ds = ds.rename_dims(lat="y", lon="x")
        # NOTE: lat/lon are 1-D in our examples, though 2-D in the example nc output

    # Store namelist settings
    ds.attrs["nml"] = str(full_config)
    ds.attrs["nml_json"] = json.dumps(full_config.todict())

    # Clean up
    if cleanup:
        shutil.rmtree(case_dir)

    return ds


def read_txt(fp: Path) -> pd.DataFrame:
    """Read canopy-app txt output file."""
    import re

    # Parse header lines
    with open(fp) as f:
        for i, line in enumerate(f):
            if i == 0:
                pattern = r" *Reference height, h\: *([0-9\.]*) m"
                m = re.match(pattern, line)
                if m is None:
                    raise ValueError(
                        f"Unexpected file format. Line {i} failed to match regex {pattern!r}."
                    )
                href = float(m.group(1))
            elif i == 1:
                pattern = r" *Number of model layers\: *([0-9]*)"
                m = re.match(pattern, line)
                if m is None:
                    raise ValueError(
                        f"Unexpected file format. Line {i} failed to match regex {pattern!r}."
                    )
                nlay = int(m.group(1))
            elif i == 2:
                # Column names (some with units)
                heads = re.split(r" {2,}", line.strip())
                names: list[str] = []
                units: dict[str, str | None] = {}
                for head in heads:
                    m = re.fullmatch(r"([a-zA-Z_ ]+) \(([^)]+)\)", head)
                    if m is None:
                        u = None
                        name = head
                    else:
                        name, u = m.groups()
                    clean_name = name.lower().replace(" ", "_")
                    names.append(clean_name)
                    units[clean_name] = u
            else:
                break
        else:
            raise ValueError(
                "Unexpected file format. Expected 3 header lines followed by data."
            )

    df = pd.read_csv(fp, index_col=False, skiprows=3, header=None, delimiter=r"\s+")
    if len(names) != len(df.columns):
        raise RuntimeError(
            f"Unexpected file format. Detected columns names {names} ({len(names)}) "
            f"are of a different number than the loaded dataframe ({len(df.columns)})."
        )
    df.columns = names
    df.attrs.update(href=href, nlay=nlay, units=units)

    return df


def sens_config(
    cases: list[dict[str, Any]],
    base_dir: Path | None = None,
    cleanup: bool = True,
    verbose: bool = False,
) -> xr.Dataset:
    """Do multiple runs with different namelist configuration,
    returning single merged dataset.
    """
    from collections import defaultdict

    if base_dir is None:
        import tempfile

        temp_dir = tempfile.TemporaryDirectory(prefix="canopy-app_")
        base_dir = Path(temp_dir.name)

    dss = []
    case_vars = defaultdict(list)
    for i, case in enumerate(cases):
        assert not ({"modlays", "modres"} & set(case["userdefs"]))

        for k, v in case["userdefs"].items():
            case_vars[k].append(v)

        if verbose:
            print(f"Running case {i+1}/{len(cases)}")
        ds = run(
            config=case,
            case_dir=base_dir / f"case_{i:0{len(str(len(cases) - 1))}}",
            cleanup=False,
            verbose=verbose,
        )
        dss.append(ds)

    ds = xr.concat(dss, dim="case", join="exact", combine_attrs="drop_conflicts")

    for k, v in case_vars.items():
        ds[k] = ("case", v)

    # Clean up
    if cleanup:
        shutil.rmtree(base_dir)

    return ds


def _k_sec(k: str) -> str:
    if k.startswith("file_"):
        sec = "filenames"
    else:
        sec = "userdefs"
    return sec


def config_cases(*, product: bool = False, **kwargs) -> list[dict[str, Any]]:
    """Helper to generate config dicts.

    Supply any namelist parameter as a keyword argument,
    with a single option or list of options
    (those with lists must be the same length, unless using `product`).
    """
    from collections import defaultdict

    if not product:
        sings = {}
        mults = {}
        for k, v in kwargs.items():
            if np.isscalar(DEFAULTS[_k_sec(k)][k]):
                sings[k] = v
            else:
                mults[k] = v

        max_len = max(len(v) for v in mults.values())
        for k, v in mults.items():
            if len(v) != max_len:
                raise ValueError(
                    "All lists should be the same length. "
                    f"Max length is {max_len} but {k!r} has length {len(v)}."
                )

        cases = []
        for i in range(max_len):
            case = defaultdict(dict)
            for k, v in sings.items():
                case[_k_sec(k)][k] = v
            for k, v in mults.items():
                case[_k_sec(k)][k] = v[i]
            cases.append(case)

    else:
        import itertools

        for k, v in kwargs.items():
            if not isinstance(v, list):
                kwargs[k] = [v]

        cases = []
        for vs in itertools.product(*kwargs.values()):
            case = defaultdict(dict)
            for k, v in zip(kwargs, vs):
                case[_k_sec(k)][k] = v
            cases.append(case)

    for i, case in enumerate(cases):
        cases[i] = dict(case)

    return cases


if __name__ == "__main__":
    ...
    # ds = run(case_dir=Path("test"), cleanup=False)
    # ds = run(
    #     config={
    #         "filenames": {"file_vars": "../input/input_variables_point.txt"},
    #         "userdefs": {"infmt_opt": 1, "nlat": 1, "nlon": 1},
    #         # "filenames": {"file_vars": "../input/gfs.t12z.20220701.sfcf000.canopy.txt"},
    #         # "userdefs": {"infmt_opt": 1},
    #     },
    #     case_dir=Path("test"),
    #     cleanup=False,
    # )

    # df = read_txt(Path("test/output/out_output_waf.txt"))

    # cases = [
    #     {
    #         "filenames": {"file_vars": "../input/input_variables_point.txt"},
    #         "userdefs": {"infmt_opt": 1, "nlat": 1, "nlon": 1, "z0ghc": 0.001},
    #     },
    #     {
    #         "filenames": {"file_vars": "../input/input_variables_point.txt"},
    #         "userdefs": {"infmt_opt": 1, "nlat": 1, "nlon": 1, "z0ghc": 0.01},
    #     },
    # ]
    # ds = sens_config(cases)

    cases = config_cases(
        file_vars="../input/input_variables_point.txt",
        infmt_opt=1,
        nlat=1,
        nlon=1,
        z0ghc=[0.001, 0.01],
        lambdars=[1.0, 1.25],
        product=True,
    )
    from pprint import pprint

    pprint(cases)
    ds = sens_config(cases)
