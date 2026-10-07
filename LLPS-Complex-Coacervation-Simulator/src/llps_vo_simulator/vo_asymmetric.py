"""Electroneutral asymmetric Voorn--Overbeek stability maps.

The asymmetric pages track polyion molecules, their monovalent counterions,
and a selectable added electrolyte.  The small-ion species are kept distinct
in the ideal-mixing entropy; all species contribute with their charge squared
to Debye--Huckel screening.  The maps report spinodals, not binodals.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi

import numpy as np


@dataclass(frozen=True)
class SaltSpec:
    key: str
    name: str
    cation_charge: int
    anion_charge: int
    cation_stoich: int
    anion_stoich: int

    @property
    def ionic_strength_factor(self) -> float:
        """I_s/c_s, where I_s = 1/2 sum_i c_i z_i^2."""
        return 0.5 * (
            self.cation_stoich * self.cation_charge**2
            + self.anion_stoich * self.anion_charge**2
        )


SALTS = {
    "NaCl": SaltSpec("NaCl", "NaCl", +1, -1, 1, 1),
    "KCl": SaltSpec("KCl", "KCl", +1, -1, 1, 1),
    "MgCl2": SaltSpec("MgCl2", "MgCl₂", +2, -1, 1, 2),
    "CaCl2": SaltSpec("CaCl2", "CaCl₂", +2, -1, 1, 2),
    "Na2SO4": SaltSpec("Na2SO4", "Na₂SO₄", +1, -2, 2, 1),
    "K2SO4": SaltSpec("K2SO4", "K₂SO₄", +1, -2, 2, 1),
    "MgSO4": SaltSpec("MgSO4", "MgSO₄", +2, -2, 1, 1),
}


def get_salt(value: str | None) -> SaltSpec:
    key = str(value or "NaCl")
    try:
        return SALTS[key]
    except KeyError as exc:
        raise ValueError(f"Unsupported added salt: {key}.") from exc


@dataclass(frozen=True)
class AsymmetricVOParameters:
    segments_p: int = 80
    segments_q: int = 80
    charge_p: float = 0.75
    charge_q: float = 0.75
    lattice_nm: float = 0.6
    bjerrum_nm: float = 0.7
    chi_polymer_water: float = 0.0

    @classmethod
    def from_mapping(cls, value: dict) -> "AsymmetricVOParameters":
        par = cls(
            segments_p=int(value.get("segments_p", 80)),
            segments_q=int(value.get("segments_q", 80)),
            charge_p=float(value.get("charge_p", 0.75)),
            charge_q=float(value.get("charge_q", 0.75)),
            lattice_nm=float(value.get("lattice_nm", 0.6)),
            bjerrum_nm=float(value.get("bjerrum_nm", 0.7)),
            chi_polymer_water=float(value.get("chi_polymer_water", 0.0)),
        )
        if not 1 <= par.segments_p <= 400 or not 1 <= par.segments_q <= 400:
            raise ValueError("Both chain lengths must be between 1 and 400 segments.")
        if not 0.05 <= par.charge_p <= 1 or not 0.05 <= par.charge_q <= 1:
            raise ValueError("Both charge-site fractions must be between 0.05 and 1.")
        if not 0.56 <= par.lattice_nm <= 1.2:
            raise ValueError("Lattice spacing must be between 0.56 and 1.2 nm.")
        if not 0.2 <= par.bjerrum_nm <= 1:
            raise ValueError("Bjerrum length must be between 0.2 and 1.0 nm.")
        if par.bjerrum_nm / par.lattice_nm > 1.8:
            raise ValueError("Bjerrum length / lattice spacing must not exceed 1.8.")
        if not -0.4 <= par.chi_polymer_water <= 0.5:
            raise ValueError("Polymer--water chi must be between -0.4 and 0.5.")
        return par

    @property
    def dh_coefficient(self) -> float:
        return (2 * np.sqrt(pi) / 3) * (self.bjerrum_nm / self.lattice_nm) ** 1.5

    @property
    def lattice_molarity(self) -> float:
        return 1 / (0.602214076 * self.lattice_nm**3)


def _xlogx(value: float) -> float:
    return value * np.log(max(value, 1e-300)) if value > 0 else 0.0


def _composition(total_polyion_m: float, salt_m: float, charge_ratio: float,
                 par: AsymmetricVOParameters, salt: SaltSpec) -> dict:
    """Build feed composition from molecule molarity and charge equivalents.

    R_q = (q_P N_P c_P)/(q_Q N_Q c_Q).  All volume fractions use the
    segment lattice molarity.  Each charged site has one monovalent
    polyion-derived counterion; added salt ions follow its formula stoichiometry.
    """
    if not np.isfinite(total_polyion_m) or not np.isfinite(salt_m):
        raise ValueError("Feed concentrations must be finite numbers.")
    if total_polyion_m < 0 or salt_m < 0:
        raise ValueError("Feed concentrations cannot be negative.")
    if not 0.1 <= charge_ratio <= 5:
        raise ValueError("Feed charge ratio Rq must be between 0.1 and 5.")

    # Resolve c_P + c_Q = C_total and the requested charge-equivalent ratio.
    p_weight = charge_ratio * par.charge_q * par.segments_q
    q_weight = par.charge_p * par.segments_p
    c_p_m = total_polyion_m * p_weight / (p_weight + q_weight)
    c_q_m = total_polyion_m - c_p_m
    p = par.segments_p * c_p_m / par.lattice_molarity
    q = par.segments_q * c_q_m / par.lattice_molarity
    x = par.charge_p * p
    y = par.charge_q * q
    salt_formula = salt_m / par.lattice_molarity
    salt_plus = salt.cation_stoich * salt_formula
    salt_minus = salt.anion_stoich * salt_formula
    water = 1 - (p + q + x + y + salt_plus + salt_minus)
    if water <= 1e-8:
        raise ValueError("This feed exceeds the lattice packing limit. Lower polymer or salt concentration.")

    ionic_screen = (
        par.charge_p * p + par.charge_q * q + x + y
        + salt.cation_charge**2 * salt_plus
        + salt.anion_charge**2 * salt_minus
    )
    charge_residual = (
        par.charge_p * p - par.charge_q * q - x + y
        + salt.cation_charge * salt_plus
        + salt.anion_charge * salt_minus
    ) * par.lattice_molarity
    cation_charge_m = y * par.lattice_molarity + salt.cation_charge * salt_plus * par.lattice_molarity
    anion_charge_m = x * par.lattice_molarity - salt.anion_charge * salt_minus * par.lattice_molarity
    return {
        "phi_p": p, "phi_q": q, "counterion_anion": x,
        "counterion_cation": y, "salt_cation": salt_plus,
        "salt_anion": salt_minus, "water": water,
        "polymer_m": total_polyion_m,
        "polymer_segment_m": (p + q) * par.lattice_molarity,
        "polycation_m": c_p_m, "polyanion_m": c_q_m,
        "salt_m": salt_m, "charge_ratio": charge_ratio,
        "counterion_anion_m": x * par.lattice_molarity,
        "counterion_cation_m": y * par.lattice_molarity,
        "salt_cation_m": salt_plus * par.lattice_molarity,
        "salt_anion_m": salt_minus * par.lattice_molarity,
        "salt_ionic_strength_m": salt.ionic_strength_factor * salt_m,
        "salt_ionic_strength_factor": salt.ionic_strength_factor,
        "polymer_positive_charge_m": par.charge_p * p * par.lattice_molarity,
        "polymer_negative_charge_m": par.charge_q * q * par.lattice_molarity,
        "mobile_cation_charge_m": cation_charge_m,
        "mobile_anion_charge_m": anion_charge_m,
        "ionic_screen_m": ionic_screen * par.lattice_molarity,
        "charge_residual_m": charge_residual,
    }


def feed_state(total_polyion_m: float, salt_m: float, charge_ratio: float,
               par: AsymmetricVOParameters, salt_type: str = "NaCl") -> dict:
    """Return a charge-neutral feed in lattice volume fractions.

    ``total_polyion_m`` is c_P + c_Q in chain-molecule molarity. ``charge_ratio``
    is the ratio of polymer charge equivalents, including chain lengths.
    ``salt_m`` is formula-unit molarity of the selected added salt.
    """
    salt = get_salt(salt_type)
    state = _composition(total_polyion_m, salt_m, charge_ratio, par, salt)
    state["salt_type"] = salt.key
    state["salt_name"] = salt.name
    state["min_eigenvalue"] = (
        spinodal_eigenvalue(total_polyion_m, salt_m, charge_ratio, par, salt.key)
        if total_polyion_m > 0 else None
    )
    if state["min_eigenvalue"] is None:
        state["local_stability"] = "not evaluated at zero polyion concentration"
    else:
        state["local_stability"] = (
            "unstable" if state["min_eigenvalue"] < 0 else "locally stable"
        )
    return state


def free_energy(phi_p: float, phi_q: float, counterion_anion: float,
                counterion_cation: float, salt_m: float,
                par: AsymmetricVOParameters, salt_type: str = "NaCl") -> float:
    """Dimensionless VO density f*a^3/(kBT) for an explicit ion inventory."""
    salt = get_salt(salt_type)
    salt_formula = salt_m / par.lattice_molarity
    salt_plus = salt.cation_stoich * salt_formula
    salt_minus = salt.anion_stoich * salt_formula
    water = 1 - (phi_p + phi_q + counterion_anion + counterion_cation + salt_plus + salt_minus)
    species = (phi_p, phi_q, counterion_anion, counterion_cation, salt_plus, salt_minus)
    if water <= 0 or any(value < 0 for value in species):
        return float("inf")
    ionic_screen = (
        par.charge_p * phi_p + par.charge_q * phi_q
        + counterion_anion + counterion_cation
        + salt.cation_charge**2 * salt_plus
        + salt.anion_charge**2 * salt_minus
    )
    mixing = (
        _xlogx(phi_p) / par.segments_p + _xlogx(phi_q) / par.segments_q
        + sum(_xlogx(value) for value in species[2:])
        + _xlogx(water)
    )
    electrostatic = -par.dh_coefficient * ionic_screen**1.5
    solvent = par.chi_polymer_water * (phi_p + phi_q) * water
    return float(mixing + electrostatic + solvent)


def spinodal_eigenvalue(total_polyion_m: float, salt_m: float, charge_ratio: float,
                       par: AsymmetricVOParameters, salt_type: str = "NaCl") -> float:
    """Smallest constrained-Hessian eigenvalue at the homogeneous feed state.

    Independent coordinates are P, Q, X-, and each present salt ion species.
    Y+ is eliminated by local electroneutrality. Eigenvalues below zero indicate
    local instability. The homogeneous map uses chain-molecule concentration.
    """
    salt = get_salt(salt_type)
    state = _composition(total_polyion_m, salt_m, charge_ratio, par, salt)
    p, q = state["phi_p"], state["phi_q"]
    x, y = state["counterion_anion"], state["counterion_cation"]
    sp, sm = state["salt_cation"], state["salt_anion"]
    water = state["water"]
    if min(p, q, x, y, water) <= 0:
        raise ValueError("Spinodal requires positive polymer, counterion and solvent fractions.")

    # Full species order: P, Q, X-, Y+, salt cation, salt anion.  Enforce
    # charge neutrality by expressing Y+ as a function of the other species.
    columns: list[tuple[int, float]] = [(0, -par.charge_p), (1, par.charge_q), (2, 1.0)]
    if sp > 0:
        columns.append((4, -float(salt.cation_charge)))
    if sm > 0:
        columns.append((5, -float(salt.anion_charge)))
    basis = np.zeros((6, len(columns)), dtype=float)
    basis[3, :3] = [-par.charge_p, par.charge_q, 1.0]
    for col, (species_index, y_gradient) in enumerate(columns[3:], start=3):
        basis[species_index, col] = 1.0
        basis[3, col] = y_gradient
    # P, Q and X- are the first three independent directions.
    basis[0, 0] = basis[1, 1] = basis[2, 2] = 1.0

    species = np.array([p, q, x, y, sp, sm], dtype=float)
    sizes = np.array([par.segments_p, par.segments_q, 1, 1, 1, 1], dtype=float)
    screen_weights = np.array([
        par.charge_p, par.charge_q, 1, 1,
        salt.cation_charge**2, salt.anion_charge**2,
    ], dtype=float)
    grad_water = -basis.sum(axis=0)
    grad_screen = screen_weights @ basis
    hessian = np.zeros((basis.shape[1], basis.shape[1]), dtype=float)
    for index, fraction in enumerate(species):
        if fraction > 0:
            gradient = basis[index]
            hessian += np.outer(gradient, gradient) / (sizes[index] * fraction)
    hessian += np.outer(grad_water, grad_water) / water

    ionic_screen = float(screen_weights @ species)
    if ionic_screen > 0:
        hessian -= (
            3 * par.dh_coefficient / (4 * np.sqrt(ionic_screen))
        ) * np.outer(grad_screen, grad_screen)
    grad_polymer = basis[0] + basis[1]
    hessian += par.chi_polymer_water * (
        np.outer(grad_polymer, grad_water) + np.outer(grad_water, grad_polymer)
    )
    return float(np.linalg.eigvalsh(hessian)[0])


def simulate_2d(par: AsymmetricVOParameters, charge_ratio: float,
                polymer_max_m: float = 0.025, salt_max_m: float = 1.0,
                polymer_count: int = 48, salt_count: int = 42,
                feed_polymer_m: float | None = None,
                feed_salt_m: float | None = None,
                salt_type: str = "NaCl") -> dict:
    salt_spec = get_salt(salt_type)
    if not 0.1 <= charge_ratio <= 5:
        raise ValueError("Feed charge ratio Rq must be between 0.1 and 5.")
    if not 0.0005 <= polymer_max_m <= 0.5 or not 0.001 <= salt_max_m <= 5:
        raise ValueError("Map limits exceed supported molecule/salt concentration range.")
    if not 8 <= polymer_count <= 120 or not 8 <= salt_count <= 120:
        raise ValueError("2D map resolution must be between 8 and 120 samples per axis.")
    polymer = np.linspace(max(polymer_max_m / polymer_count * 0.05, 1e-8), polymer_max_m, polymer_count)
    salt_values = np.linspace(0.0, salt_max_m, salt_count)
    eigen = np.full((salt_count, polymer_count), np.nan)
    for iy, salt_m in enumerate(salt_values):
        for ix, polymer_m in enumerate(polymer):
            try:
                eigen[iy, ix] = spinodal_eigenvalue(
                    polymer_m, salt_m, charge_ratio, par, salt_spec.key
                )
            except ValueError:
                pass
    selected_feed = None
    if feed_polymer_m is not None and feed_salt_m is not None:
        selected_feed = feed_state(
            float(feed_polymer_m), float(feed_salt_m), charge_ratio, par, salt_spec.key
        )
    return {
        "parameters": par.__dict__, "charge_ratio": charge_ratio,
        "polymer_m": polymer.tolist(), "salt_m": salt_values.tolist(),
        "min_eigenvalue": [[None if not np.isfinite(v) else float(v) for v in row] for row in eigen],
        "unstable": [[bool(np.isfinite(v) and v < 0) for v in row] for row in eigen],
        "lattice_molarity": par.lattice_molarity,
        "dh_coefficient": par.dh_coefficient,
        "salt_type": salt_spec.key, "salt_name": salt_spec.name,
        "salt_ionic_strength_factor": salt_spec.ionic_strength_factor,
        "feed": selected_feed,
        "note": (
            "Spinodal: the smallest constrained-Hessian eigenvalue crosses zero. "
            "Polyion concentration is c_P+c_Q (chain molecules); salt concentration "
            "is formula-unit molarity. This local-stability boundary is not the binodal."
        ),
    }


def simulate_high_dim(par: AsymmetricVOParameters, rq_min: float = 0.25,
                      rq_max: float = 4.0, polymer_max_m: float = 0.025,
                      salt_max_m: float = 1.0, ratio_count: int = 19,
                      polymer_count: int = 24, salt_count: int = 61,
                      salt_type: str = "NaCl", polymer_min_m: float = 0.0001,
                      salt_min_m: float = 0.0) -> dict:
    salt_spec = get_salt(salt_type)
    if not 0.1 <= rq_min < rq_max <= 5:
        raise ValueError("Scan limits must satisfy 0.1 <= Rq min < Rq max <= 5.")
    if not 0 < polymer_min_m < polymer_max_m <= 0.5 or not 0 <= salt_min_m < salt_max_m <= 5:
        raise ValueError("Map limits exceed supported molecule/salt concentration range.")
    if not 7 <= ratio_count <= 51 or not 8 <= polymer_count <= 100 or not 10 <= salt_count <= 200:
        raise ValueError("3D resolution is outside supported bounds.")
    ratios = np.linspace(rq_min, rq_max, ratio_count)
    polymers = np.linspace(polymer_min_m, polymer_max_m, polymer_count)
    salts = np.linspace(salt_min_m, salt_max_m, salt_count)
    boundary: list[list[float | None]] = []
    volume: list[dict] = []
    for ratio in ratios:
        ratio_row = []
        for polymer_m in polymers:
            vals = []
            for salt_m in salts:
                try:
                    vals.append(spinodal_eigenvalue(
                        polymer_m, salt_m, float(ratio), par, salt_spec.key
                    ))
                except ValueError:
                    vals.append(float("nan"))
            crossings = []
            for i in range(len(salts) - 1):
                left, right = vals[i], vals[i + 1]
                if np.isfinite(left) and np.isfinite(right) and left < 0 <= right:
                    frac = -left / (right - left) if right != left else 0.0
                    crossings.append(float(salts[i] + frac * (salts[i + 1] - salts[i])))
            threshold = crossings[0] if crossings else None
            ratio_row.append(threshold)
            if threshold is not None:
                volume.append({
                    "charge_ratio": float(ratio),
                    "polymer_m": float(polymer_m),
                    "salt_m": threshold,
                    "salt_ionic_strength_m": salt_spec.ionic_strength_factor * threshold,
                })
        boundary.append(ratio_row)
    return {
        "parameters": par.__dict__, "charge_ratio": ratios.tolist(),
        "polymer_m": polymers.tolist(), "salt_boundary_m": boundary,
        "points": volume, "salt_min_m": salt_min_m, "salt_max_m": salt_max_m,
        "salt_type": salt_spec.key, "salt_name": salt_spec.name,
        "salt_ionic_strength_factor": salt_spec.ionic_strength_factor,
        "salt_ionic_strength_boundary_m": [
            [None if value is None else salt_spec.ionic_strength_factor * value for value in row]
            for row in boundary
        ],
        "axes": {
            "x": "polymer charge-equivalent ratio R_q = q_P N_P c_P / (q_Q N_Q c_Q)",
            "y": "total polyion molecule concentration c_P + c_Q (M)",
            "z": f"{salt_spec.name} formula-unit concentration c_s (M)",
        },
        "note": (
            "Three-dimensional spinodal surface: polymer charge-equivalent ratio, "
            "total polyion molecule concentration, and selected added-salt formula "
            "concentration. Salt ionic strength is I_s = "
            f"{salt_spec.ionic_strength_factor:g} c_s."
        ),
    }
