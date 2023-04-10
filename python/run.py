"""
Run canopy-app using the namelist and load nc data.
"""
from __future__ import annotations

from pathlib import Path

import f90nml

HERE = Path(__file__).parent
assert HERE.name == "python"
assert HERE.parent.name == "canopy-app"
REPO = HERE.parent

with open(REPO / "input" / "namelist.canopy") as f:
    DEFAULTS = f90nml.read(f)

print(DEFAULTS)


def run(config: dict):
    ...
