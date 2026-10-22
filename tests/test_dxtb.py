import pytest
import torch

from config import ValidatorConfig
from dxtb_backend import DxtbError, DxtbRunner
from gfn2 import Gfn2Runner
from validator import RouteValidator

# GFN1 needs no libcint, so it runs here. GFN2 in dxtb needs libcint, and its test runs only when that is present.
requires_gfn1 = pytest.mark.skipif(not DxtbRunner("GFN1").available(), reason="dxtb GFN1 is not available")


@requires_gfn1
def test_relaxation_lowers_the_energy_below_the_starting_geometry():
    runner = DxtbRunner("GFN1")
    numbers, positions, charge, _ = Gfn2Runner.geometry("CC")
    start = runner._calculator(numbers).get_energy(torch.tensor(positions, dtype=torch.double), chrg=charge)
    assert runner.energy("CC") <= float(start) + 1e-8


@requires_gfn1
def test_h2_enthalpy_and_free_energy_are_consistent_with_entropy():
    energy, enthalpy, free = DxtbRunner("GFN1").thermo("[H][H]")
    # H - E is ZPE plus thermal terms, roughly 30 to 50 kJ/mol for H2 in GFN1
    assert 30 < (enthalpy - energy) * 2625.5 < 50
    # G = H - TS, and TS is about 40 kJ/mol for H2 at 298 K, so G sits below H
    assert free < enthalpy


@requires_gfn1
def test_water_has_no_imaginary_modes():
    # thermochemistry raises an EnergyError if any mode is imaginary
    energy, enthalpy, free = DxtbRunner("GFN1").thermo("O")
    assert free < enthalpy


@requires_gfn1
def test_dxtb_validator_reports_dE_dH_dG_for_a_hydrogenation_step(tmp_path):
    validator = RouteValidator(ValidatorConfig(backend="dxtb", dxtb_method="GFN1", workdir=str(tmp_path)))
    detail = validator.step_detail("CC", {"C=C": 1, "[H][H]": 1})
    assert set(detail) == {"dE", "dH", "dG"}
    assert detail["dE"] < 0
    assert detail["dG"] < 0  # hydrogenation stays favourable after the thermal terms
    assert abs(detail["dH"] - detail["dE"]) > 1  # the thermal terms are applied


def test_gfn2_in_dxtb_reports_missing_libcint_clearly():
    runner = DxtbRunner("GFN2")
    if runner.available():
        pytest.skip("libcint is installed here, so GFN2 runs")
    with pytest.raises(DxtbError, match="libcint"):
        runner.energy("[H][H]")


def test_unknown_dxtb_method_is_rejected():
    with pytest.raises(ValueError, match="unknown dxtb method"):
        DxtbRunner("GFN0")


def test_dxtb_backend_allows_free_energy_criterion(tmp_path):
    validator = RouteValidator(ValidatorConfig(backend="dxtb", criterion="dg", workdir=str(tmp_path)))
    assert validator.criterion_label == "dG"

