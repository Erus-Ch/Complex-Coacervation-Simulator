"""Symmetric Voorn--Overbeek coexistence calculation, using the symmetric three-component reduction.

The conserved coordinates are the total polymer segment volume fraction p and
the volume fraction n of *each* mobile ion sign. The two polyions have equal
segment fractions p/2, equal N and charge-site fraction q. Their charges cancel
locally; n includes counterions and added 1:1 salt ions. This is the symmetric
three-component reduction discussed by Sing and Perry (Soft Matter 2020).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi

import numpy as np
from scipy.optimize import root
from scipy.special import expit


@dataclass(frozen=True)
class VOParameters:
    segments: int = 80
    charge_fraction: float = 0.75
    lattice_nm: float = 0.6
    bjerrum_nm: float = 0.7
    chi_polymer_water: float = 0.0

    @classmethod
    def from_mapping(cls, value: dict) -> "VOParameters":
        result = cls(
            segments=int(value.get("segments", 80)),
            charge_fraction=float(value.get("charge_fraction", 0.75)),
            lattice_nm=float(value.get("lattice_nm", 0.6)),
            bjerrum_nm=float(value.get("bjerrum_nm", 0.7)),
            chi_polymer_water=float(value.get("chi_polymer_water", 0.0)),
        )
        if not 1 <= result.segments <= 400:
            raise ValueError("N must be between 1 and 400 lattice segments.")
        if not 0.1 <= result.charge_fraction <= 1:
            raise ValueError("Charge-site fraction must be between 0.1 and 1.")
        if not 0.56 <= result.lattice_nm <= 1.2:
            raise ValueError("Lattice spacing must be between 0.56 and 1.2 nm.")
        if not 0.2 <= result.bjerrum_nm <= 1.0:
            raise ValueError("Bjerrum length must be between 0.2 and 1.0 nm.")
        if result.bjerrum_nm / result.lattice_nm > 1.8:
            raise ValueError("Bjerrum length / lattice spacing must not exceed 1.8 in this solver.")
        if not -0.4 <= result.chi_polymer_water <= 0.5:
            raise ValueError("Polymer--water chi must be between -0.4 and 0.5.")
        return result

    @property
    def dh_coefficient(self) -> float:
        # v0*kappa^3/(12*pi); kappa^2=4*pi*l_B*(q*p+2*n)/v0.
        return (2 * np.sqrt(pi) / 3) * (self.bjerrum_nm / self.lattice_nm) ** 1.5

    @property
    def lattice_molarity(self) -> float:
        # One lattice site per a^3 nm^3, converted to mol/L.
        return 1 / (0.602214076 * self.lattice_nm**3)


def _xlogx(x: np.ndarray) -> np.ndarray:
    return np.where(x > 0, x * np.log(np.maximum(x, 1e-300)), 0.0)


def free_energy(p: np.ndarray, n: np.ndarray, par: VOParameters) -> np.ndarray:
    """Reduced Helmholtz free-energy density f*v0/(k_B*T)."""
    p, n = np.asarray(p), np.asarray(n)
    w = 1 - p - 2 * n
    ionic_sites = par.charge_fraction * p + 2 * n
    if np.any(p < 0) or np.any(n < 0) or np.any(w <= 0):
        raise ValueError("Volume fractions must satisfy p>=0, n>=0, p+2n<1.")
    mixing = _xlogx(p / 2) * 2 / par.segments + 2 * _xlogx(n) + _xlogx(w)
    dh = -par.dh_coefficient * ionic_sites**1.5
    return mixing + dh + par.chi_polymer_water * p * w


def chemical_potentials(p: np.ndarray, n: np.ndarray, par: VOParameters) -> tuple[np.ndarray, np.ndarray]:
    p, n = np.asarray(p), np.asarray(n)
    w = 1 - p - 2 * n
    root_ionic = np.sqrt(np.maximum(par.charge_fraction * p + 2 * n, 1e-300))
    mu_p = (
        (np.log(np.maximum(p / 2, 1e-300)) + 1) / par.segments
        - np.log(np.maximum(w, 1e-300)) - 1
        - 1.5 * par.dh_coefficient * par.charge_fraction * root_ionic
        + par.chi_polymer_water * (w - p)
    )
    mu_n = (
        2 * np.log(np.maximum(n / w, 1e-300))
        - 3 * par.dh_coefficient * root_ionic
        - 2 * par.chi_polymer_water * p
    )
    return mu_p, mu_n


def _ion_fraction_at_mu(p: np.ndarray, mu_n: float, par: VOParameters) -> np.ndarray:
    """Solve the salt chemical-potential equation on the physical domain.

    The selected parameter bounds keep the n direction convex for the plotted
    states. Bisection is vectorized so a complete binodal can update promptly.
    """
    p = np.asarray(p, dtype=float)
    lo = np.full_like(p, 1e-13)
    hi = (1 - p) / 2 - 1e-13
    for _ in range(45):
        mid = (lo + hi) / 2
        trial = chemical_potentials(p, mid, par)[1]
        lo = np.where(trial < mu_n, mid, lo)
        hi = np.where(trial >= mu_n, mid, hi)
    return (lo + hi) / 2


def _lower_hull_gap(p: np.ndarray, g: np.ndarray) -> tuple[int, int] | None:
    hull: list[int] = []
    for k in range(len(p)):
        while len(hull) >= 2:
            i, j = hull[-2], hull[-1]
            slope_before = (g[j] - g[i]) / (p[j] - p[i])
            slope_after = (g[k] - g[j]) / (p[k] - p[j])
            if slope_before < slope_after:
                break
            hull.pop()
        hull.append(k)
    candidate = max(zip(hull[:-1], hull[1:]), key=lambda pair: pair[1] - pair[0])
    i, j = candidate
    if j - i < 3:
        return None
    line = g[i] + (g[j] - g[i]) * (p[i:j + 1] - p[i]) / (p[j] - p[i])
    if np.max(g[i:j + 1] - line) < 1e-10:
        return None
    return candidate


def _binodal_at_mu(mu_n: float, p_grid: np.ndarray, par: VOParameters) -> dict | None:
    n = _ion_fraction_at_mu(p_grid, mu_n, par)
    g = free_energy(p_grid, n, par) - mu_n * n
    gap = _lower_hull_gap(p_grid, g)
    if gap is None:
        return None
    i, j = gap
    p_a, p_b = float(p_grid[i]), float(p_grid[j])
    n_a, n_b = float(n[i]), float(n[j])

    # The hull locates the two branches; refine their compositions until the
    # polymer chemical potentials and semi-grand tangent intercepts coincide.
    # Log/logit coordinates allow a polymer-poor branch very close to zero.
    def tangent_residual(coordinates: np.ndarray) -> np.ndarray:
        a = float(np.exp(np.clip(coordinates[0], -700, -1e-7)))
        b = float(expit(coordinates[1]))
        if a >= b or b >= 0.999:
            return np.array([1 + a - b, 1 + a - b])
        polymer = np.array([a, b])
        ions = _ion_fraction_at_mu(polymer, mu_n, par)
        semigrand = free_energy(polymer, ions, par) - mu_n * ions
        mu_polymer = chemical_potentials(polymer, ions, par)[0]
        intercept = semigrand - mu_polymer * polymer
        return np.array([mu_polymer[0] - mu_polymer[1], intercept[0] - intercept[1]])

    initial = np.array([np.log(p_a), np.log(p_b / (1 - p_b))])
    solution = root(tangent_residual, initial, tol=1e-10)
    if solution.success and np.max(np.abs(tangent_residual(solution.x))) < 1e-7:
        candidate_a = float(np.exp(np.clip(solution.x[0], -700, -1e-7)))
        candidate_b = float(expit(solution.x[1]))
        if candidate_b - candidate_a > 0.1 * (p_b - p_a) and candidate_b < 0.94:
            p_a, p_b = candidate_a, candidate_b
            n_a, n_b = map(float, _ion_fraction_at_mu(np.array([p_a, p_b]), mu_n, par))
    return {
        "dilute": [p_a, n_a],
        "dense": [p_b, n_b],
        "mu_ion_pair": float(mu_n),
        "gap": p_b - p_a,
        "grand_excess_max": float(np.max(g[i:j + 1] - (g[i] + (g[j] - g[i]) * (p_grid[i:j + 1] - p_grid[i]) / (p_grid[j] - p_grid[i])))),
    }


def _profile(pair: dict, par: VOParameters) -> dict:
    a, b = pair["dilute"], pair["dense"]
    upper = min(0.96, b[0] * 1.08)
    knee = min(max(a[0] * 10, 1e-5), upper * 0.5)
    p = np.unique(np.r_[
        np.geomspace(max(a[0] * 0.3, 1e-20), knee, 80),
        np.linspace(knee, upper, 300), a[0], b[0],
    ])
    n = _ion_fraction_at_mu(p, pair["mu_ion_pair"], par)
    g = free_energy(p, n, par) - pair["mu_ion_pair"] * n
    ga = float(free_energy(np.array(a[0]), np.array(a[1]), par) - pair["mu_ion_pair"] * a[1])
    gb = float(free_energy(np.array(b[0]), np.array(b[1]), par) - pair["mu_ion_pair"] * b[1])
    tangent = ga + (gb - ga) * (p - a[0]) / (b[0] - a[0])
    return {
        "polymer_m": [float(f"{value:.12g}") for value in p * par.lattice_molarity],
        "grand_above_tangent": [float(f"{value:.12g}") for value in g - tangent],
        "dilute_polymer_m": a[0] * par.lattice_molarity,
        "dense_polymer_m": b[0] * par.lattice_molarity,
        "ion_activity_mu": pair["mu_ion_pair"],
    }


def simulate(par: VOParameters) -> dict:
    """Trace the binodal by convexifying the salt semi-grand free energy."""
    p_grid = np.unique(np.r_[np.geomspace(1e-7, 0.012, 92), np.linspace(0.0125, 0.94, 440)])
    reservoir_n = np.r_[np.geomspace(1e-6, 0.01, 26), np.linspace(0.011, 0.43, 105)]
    pairs = []
    last_present_n0 = None
    first_absent_after_present = None
    for n0 in reservoir_n:
        mu = float(chemical_potentials(np.array(1e-10), np.array(n0), par)[1])
        pair = _binodal_at_mu(mu, p_grid, par)
        if pair:
            pairs.append(pair)
            last_present_n0 = n0
        elif last_present_n0 is not None and first_absent_after_present is None:
            first_absent_after_present = n0
    if last_present_n0 is not None and first_absent_after_present is not None:
        lo, hi = last_present_n0, first_absent_after_present
        for _ in range(15):
            mid = (lo + hi) / 2
            mu = float(chemical_potentials(np.array(1e-10), np.array(mid), par)[1])
            pair = _binodal_at_mu(mu, p_grid, par)
            if pair:
                lo = mid
                pairs.append(pair)
            else:
                hi = mid
    if not pairs:
        return {
            "parameters": par.__dict__,
            "dh_coefficient": par.dh_coefficient,
            "lattice_molarity": par.lattice_molarity,
            "pairs": [],
            "profile": None,
            "message": "No two-phase region was found in the scanned VO composition range for these parameters.",
        }
    # The hull can pick up a boundary artifact if the dense branch reaches the
    # p-grid limit. Such points are excluded rather than drawn as coexistence.
    pairs = [pair for pair in pairs if pair["dense"][0] < 0.935]
    if not pairs:
        return {
            "parameters": par.__dict__, "dh_coefficient": par.dh_coefficient,
            "lattice_molarity": par.lattice_molarity, "pairs": [], "profile": None,
            "message": "The dense phase lies beyond the model's plotted composition domain. Reduce the electrostatic strength.",
        }
    c0 = par.lattice_molarity
    rendered = [
        {
            "dilute": [float(f"{pair['dilute'][0] * c0:.12g}"), float(f"{pair['dilute'][1] * c0:.12g}")],
            "dense": [float(f"{pair['dense'][0] * c0:.12g}"), float(f"{pair['dense'][1] * c0:.12g}")],
            "mu_ion_pair": pair["mu_ion_pair"],
            "grand_excess_max": pair["grand_excess_max"],
        }
        for pair in pairs
    ]
    representative = pairs[len(pairs) // 2]
    return {
        "parameters": par.__dict__,
        "dh_coefficient": par.dh_coefficient,
        "lattice_molarity": c0,
        "pairs": rendered,
        "profile": _profile(representative, par),
        "critical_estimate": [round((pairs[-1]["dilute"][0] + pairs[-1]["dense"][0]) * c0 / 2, 6), round((pairs[-1]["dilute"][1] + pairs[-1]["dense"][1]) * c0 / 2, 6)],
        "salt_partition_sign": "coacervate" if np.mean([x["dense"][1] - x["dilute"][1] for x in pairs]) > 0 else "supernatant",
        "message": "VO binodal calculated from a common tangent of the salt semi-grand free energy.",
    }


def profile_at_mu(par: VOParameters, mu_n: float) -> dict:
    p_grid = np.unique(np.r_[np.geomspace(1e-7, 0.012, 92), np.linspace(0.0125, 0.94, 440)])
    pair = _binodal_at_mu(mu_n, p_grid, par)
    if pair is None:
        raise ValueError("No coexistence at this ion chemical potential.")
    return _profile(pair, par)
