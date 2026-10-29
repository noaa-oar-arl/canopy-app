"""Compile and exercise regional MEGAN factors using only Python and gfortran."""

import itertools
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPECIES = (
    "ISOP",
    "MYRC",
    "SABI",
    "LIMO",
    "CARE",
    "OCIM",
    "BPIN",
    "APIN",
    "MONO",
    "FARN",
    "CARY",
    "SESQ",
    "MBOL",
    "METH",
    "ACET",
    "CO",
    "BVOC",
    "SVOC",
    "OVOC",
)
LATITUDES = (0.0,) + tuple(
    sign * latitude
    for sign in (-1, 1)
    for latitude in (
        10.0,
        23.499999,
        23.5,
        23.500001,
        45.0,
        66.499999,
        66.5,
        66.500001,
        80.0,
        90.0,
    )
)
VEGETATION = tuple(range(21)) + (99,)
VALID_VEGETATION = tuple(range(1, 13)) + (14, 18, 19)

DRIVER = """
program regional_factors
    use canopy_const_mod, only: rk
    use canopy_bioparm_mod, only: canopy_biop
    use canopy_bioemi_mod, only: canopy_bio
    implicit none
    integer :: species, lu, vtype, status
    real(rk) :: lat, p(23), emission(2)

    do
        read(*, *, iostat=status) species, lu, vtype, lat
        if (status < 0) exit
        if (status /= 0) stop 1
        call canopy_biop(species, lu, vtype, lat, &
            p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), &
            p(10), p(11), p(12), p(13), p(14), p(15), p(16), p(17), &
            p(18), p(19), p(20), p(21), p(22), p(23))
        select case (vtype)
        case (1:12, 14, 18:19)
            call canopy_bio(zk=[1.0_rk, 2.0_rk], fclai=[0.5_rk, 1.0_rk], &
                fch=2.0_rk, lai=2.0_rk, fsun=[0.5_rk, 0.5_rk], &
                ppfd_sun=[500.0_rk, 500.0_rk], ppfd_shade=[100.0_rk, 100.0_rk], &
                tleaf_sun=[303.0_rk, 303.0_rk], tleaf_shade=[300.0_rk, 300.0_rk], &
                ppfd24_sun=[500.0_rk, 500.0_rk], ppfd24_shade=[100.0_rk, 100.0_rk], &
                tleaf24_ave=[300.0_rk, 300.0_rk], &
                ppfd240_sun=[500.0_rk, 500.0_rk], ppfd240_shade=[100.0_rk, 100.0_rk], &
                tleaf240_ave=[300.0_rk, 300.0_rk], tka=[300.0_rk, 300.0_rk], &
                dswrf=500.0_rk, temp2=300.0_rk, lu_opt=lu, vtype=vtype, lat=lat, &
                modres=1.0_rk, cce=0.21_rk, vert=0, co2opt=1, co2set=400.0_rk, &
                leafageopt=1, pastlai=2.0_rk, currentlai=2.0_rk, tsteplai=1.0_rk, &
                lossopt=0, lossset=1.0_rk, lossind=0, lifetime=3600.0_rk, ustar=0.5_rk, &
                soimopt=1, soim1=0.3_rk, soim2=0.3_rk, soim3=0.3_rk, soim4=0.3_rk, &
                soid1=5.0_rk, soid2=25.0_rk, soid3=70.0_rk, soid4=150.0_rk, &
                wilt=0.1_rk, aqopt=2, w126_set=30.0_rk, w126_ref=30.0_rk, &
                htopt=1, ltopt=1, hwopt=1, daily_maxt2=303.0_rk, daily_mint2=290.0_rk, &
                daily_maxws10=5.0_rk, modlays=2, emi_ind=species, emi_out=emission)
            write(*, '(*(ES25.16E3,1X))') p, emission
        case default
            write(*, '(ES25.16E3)') p(1)
        end select
    end do
end program regional_factors
"""


class RegionalFactorsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory(prefix="canopy-latitude-")
        cls.addClassCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        driver = directory / "regional_factors.F90"
        driver.write_text(DRIVER, encoding="ascii")
        cls.executable = directory / "regional_factors"
        sources = (
            "canopy_const_mod.F90",
            "canopy_utils_mod.F90",
            "canopy_tleaf_mod.F90",
            "canopy_bioparm_mod.F90",
            "canopy_bioemi_mod.F90",
        )
        subprocess.run(
            [
                "gfortran",
                "-ffree-line-length-none",
                "-fcheck=all",
                "-ffpe-trap=invalid,zero,overflow",
                *(str(ROOT / "src" / source) for source in sources),
                str(driver),
                "-o",
                str(cls.executable),
            ],
            cwd=directory,
            check=True,
            capture_output=True,
            text=True,
        )
        cls.cases = tuple(itertools.product(range(1, 20), (0, 1), VEGETATION, LATITUDES))
        result = subprocess.run(
            [str(cls.executable)],
            input="".join(f"{s} {lu} {v} {lat}\n" for s, lu, v, lat in cls.cases),
            check=True,
            capture_output=True,
            text=True,
        )
        rows = result.stdout.splitlines()
        if len(rows) != len(cls.cases):
            raise AssertionError(f"Expected {len(cls.cases)} results, got {len(rows)}")
        cls.results = dict(
            zip(cls.cases, (tuple(map(float, row.split())) for row in rows))
        )
        # Read the scientific constants, not the implementation's mapping expressions.
        source = (ROOT / "src" / "canopy_bioparm_mod.F90").read_text()
        cls.factors = {
            (species, int(index)): float(value)
            for index, species, value in re.findall(
                r"PARAMETER\s*::\s*EF(\d+)_(\w+)\s*=\s*([\d.]+)_rk", source
            )
        }

    @staticmethod
    def factor_indices(vtype, lat):
        tropical, boreal = abs(lat) <= 23.5, abs(lat) >= 66.5
        return {
            1: (2,) if boreal else (1,),
            2: (4,) if tropical else (5,),
            3: (3,),
            4: (6,) if tropical else (8,) if boreal else (7,),
            5: (4, 6) if tropical else (2, 3, 8) if boreal else (1, 5, 7),
            **dict.fromkeys((6, 7, 18, 19), (11,) if boreal else (9, 10)),
            **dict.fromkeys((8, 9, 10, 11), (12, 13, 14)),
            12: (15,),
            14: (15,),
        }.get(vtype, ())

    def test_all_species_regions_and_vegetation_types(self):
        self.assertEqual(len(self.factors), 19 * 15)
        for (species, lu, vtype, lat), values in self.results.items():
            with self.subTest(species=species, lu=lu, vtype=vtype, lat=lat):
                indices = self.factor_indices(vtype, lat)
                expected = (
                    sum(self.factors[SPECIES[species - 1], i] for i in indices)
                    / len(indices)
                    if indices
                    else 0.0
                )
                self.assertAlmostEqual(values[0], expected, places=10)
                self.assertEqual(len(values), 25 if indices else 1)

    def test_other_parameters_remain_latitude_independent(self):
        for species, lu, vtype in itertools.product(
            range(1, 20), (0, 1), VALID_VEGETATION
        ):
            reference = self.results[species, 0, vtype, 0.0][1:23]
            for lat in LATITUDES:
                with self.subTest(species=species, lu=lu, vtype=vtype, lat=lat):
                    self.assertEqual(
                        self.results[species, lu, vtype, lat][1:23], reference
                    )

    def test_canopy_emissions_follow_regional_factors(self):
        for species, lu, vtype in itertools.product(
            range(1, 20), (0, 1), VALID_VEGETATION
        ):
            reference = self.results[species, lu, vtype, 0.0]
            for lat in LATITUDES:
                values = self.results[species, lu, vtype, lat]
                with self.subTest(species=species, lu=lu, vtype=vtype, lat=lat):
                    for layer in (23, 24):
                        self.assertGreater(reference[layer], 0.0)
                        self.assertAlmostEqual(
                            values[layer] / reference[layer],
                            values[0] / reference[0],
                            places=10,
                        )

    def test_invalid_latitudes_report_error(self):
        for lu, lat in itertools.product(
            (0, 1), ("-90.01", "90.01", "NaN", "Inf", "-Inf")
        ):
            with self.subTest(lu=lu, lat=lat):
                result = subprocess.run(
                    [str(self.executable)],
                    input=f"1 {lu} 1 {lat}\n",
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 2)
                self.assertIn(
                    "Invalid latitude for biogenic emission factors", result.stdout
                )


if __name__ == "__main__":
    unittest.main()
