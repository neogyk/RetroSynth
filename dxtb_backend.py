"""dxtb backend: GFN xTB energies, gradients and Hessians in PyTorch (the dxtb package).

Relaxation uses dxtb's autograd gradients with scipy L-BFGS. The Hessian is dxtb's own
get_hessian. Enthalpy and free energy use the same RRHO code as the GFN2-tblite backend.

GFN2 in dxtb needs the libcint integral library. Without it DxtbRunner("GFN2") raises DxtbError
when it starts. GFN1 needs no libcint, so the tests use it to check the plumbing.
"""

import numpy as np
from rdkit import Chem
from scipy.optimize import minimize

from gfn2 import Gfn2Runner
from thermo import (
    HARTREE_J, IMAGINARY_CM, NA, inertia_eigenvalues, is_linear, mass_weighted_frequencies, rrho_thermo,
)

METHODS = ("GFN1", "GFN2")


def _calculator_class(method):
    # Imported here, not at module load: torch and tblite each bring an OpenMP runtime,
    # and the default GFN2-tblite pipeline must not load torch.
    import dxtb.calculators as calculators

    return getattr(calculators, f"{method}Calculator")


class DxtbError(RuntimeError):
    pass


class DxtbRunner:
    """Relaxes molecules with dxtb and returns electronic energies, enthalpies and free energies."""

    def __init__(self, method="GFN2"):
        if method not in METHODS:
            raise ValueError(f"unknown dxtb method {method!r}; use one of {list(METHODS)}")
        self.method = method
        self._relaxed = {}
        self._thermo = {}

    def available(self):
        try:
            self._calculator(np.array([1, 1]))
        except DxtbError:
            return False
        return True

    def _calculator(self, numbers):
        import torch

        try:
            return _calculator_class(self.method)(torch.as_tensor(numbers), dtype=torch.double)
        except Exception as exc:  # dxtb raises RuntimeError when a required library is missing
            raise DxtbError(f"dxtb {self.method} could not start: {exc}") from exc

    def relax(self, smiles):
        """Relaxed structure: (numbers, positions in Bohr, charge, unpaired, energy). Cached."""
        key = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        if key in self._relaxed:
            return self._relaxed[key]
        numbers, positions, charge, unpaired = Gfn2Runner.geometry(key)
        import torch

        calc = self._calculator(numbers)
        shape = positions.shape
        spin = unpaired or None

        def objective(x):
            p = torch.tensor(x.reshape(shape), dtype=torch.double, requires_grad=True)
            energy = calc.get_energy(p, chrg=charge, spin=spin)
            grad, = torch.autograd.grad(energy, p)
            return float(energy.detach()), grad.numpy().ravel()

        result = minimize(
            objective, positions.ravel(), jac=True, method="L-BFGS-B",
            options={"gtol": 1e-4, "maxiter": 500},
        )
        if not result.success:
            raise DxtbError(f"dxtb geometry optimisation did not converge: {result.message}")
        self._relaxed[key] = (numbers, result.x.reshape(shape), charge, unpaired, float(result.fun))
        return self._relaxed[key]

    def energy(self, smiles):
        """Electronic energy (Hartree) of the relaxed structure."""
        return self.relax(smiles)[4]

    def hessian(self, numbers, positions, charge, unpaired):
        """Cartesian Hessian (Hartree/Bohr^2) from dxtb."""
        import torch

        calc = self._calculator(numbers)
        p = torch.tensor(positions, dtype=torch.double, requires_grad=True)
        hess = calc.get_hessian(p, chrg=charge, spin=unpaired or None)
        n = positions.size
        return hess.detach().reshape(n, n).numpy()

    def thermo(self, smiles):
        """(electronic energy, enthalpy H, free energy G) in Hartree at 298.15 K and 1 atm. Cached."""
        key = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        if key not in self._thermo:
            numbers, positions, charge, unpaired, energy = self.relax(key)
            hess = self.hessian(numbers, positions, charge, unpaired)
            masses = np.array([Chem.GetPeriodicTable().GetAtomicWeight(int(z)) for z in numbers])
            positions3 = positions.reshape(-1, 3)
            linear = is_linear(inertia_eigenvalues(positions3, masses))
            freqs = mass_weighted_frequencies(hess, masses, linear)
            thermo = rrho_thermo(positions3, masses, freqs)
            if thermo.n_imaginary:
                raise DxtbError(
                    f"{key} has an imaginary mode below {IMAGINARY_CM:.0f} cm^-1; the structure is not a minimum"
                )
            enthalpy = energy + thermo.h_corr / (NA * HARTREE_J)
            free = energy + thermo.g_corr / (NA * HARTREE_J)
            self._thermo[key] = (energy, enthalpy, free)
        return self._thermo[key]

    def enthalpy(self, smiles):
        return self.thermo(smiles)[1]

    def free_energy(self, smiles):
        return self.thermo(smiles)[2]
