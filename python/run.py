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
import xarray as xr

HERE = Path(__file__).parent.absolute()
assert HERE.name == "python"
assert HERE.parent.name == "canopy-app"
REPO = HERE.parent

with open(REPO / "input" / "namelist.canopy") as f:
    DEFAULTS = f90nml.read(f)

# Make input paths absolute in default config
for k, v in DEFAULTS["filenames"].items():
    p0 = Path(v)
    if not p0.is_absolute():
        p = REPO / p0
    else:
        p = p0
    DEFAULTS["filenames"][k] = p.as_posix()


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


def run(
    *,
    config: dict[str, Any] | None = None,
    case_dir: Path | None = None,
    cleanup: bool = True,
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
    full_config: f90nml.Namelist
    if config is None:
        full_config = DEFAULTS.copy()
    else:
        full_config = f90nml.Namelist({**DEFAULTS, **config})
    output_dir = case_dir / "output"
    output_dir.mkdir(exist_ok=True)
    ofp_stem = output_dir / "out"
    full_config["filenames"]["file_out"] = ofp_stem.relative_to(case_dir).as_posix()
    print(full_config)

    # Write namelist
    config_dir = case_dir / "input"
    config_dir.mkdir(exist_ok=True)
    full_config.write(config_dir / "namelist.canopy", force=True)

    # Run canopy-app
    exe = REPO / "canopy"
    if not exe.is_file():
        raise RuntimeError("compile canopy-app first")
    with out_and_back(case_dir):
        subprocess.run([exe])

    # Load nc
    ds = xr.open_dataset(ofp_stem.with_suffix(".nc"))

    # Store namelist settings
    ds.attrs["nml"] = str(full_config)
    ds.attrs["nml_json"] = json.dumps(full_config.todict())

    # Clean up
    if cleanup:
        shutil.rmtree(case_dir)

    return ds


if __name__ == "__main__":
    # run()
    ds = run(case_dir=Path("test"), cleanup=False)
