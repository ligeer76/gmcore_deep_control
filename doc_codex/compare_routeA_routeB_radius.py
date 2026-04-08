#!/usr/bin/env python3
"""Compare Route A and Route B radius inversion formulas.

The comparison uses a common thermodynamic integral

    X = - \int alpha d(pressure-like variable)

and defines an equivalent height h = X / G so that both routes are compared
at the same cumulative integral strength. This avoids the trivial result that
each route recovers r=a+z exactly when evaluated on its own self-consistent
vertical coordinate.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def route_b_radius(h_eq_m: np.ndarray, a_m: float) -> np.ndarray:
    """Wood-style cubic inversion for Route B."""
    return (a_m**3 + 3.0 * a_m**2 * h_eq_m) ** (1.0 / 3.0)


def route_a_radius(h_eq_m: np.ndarray, a_m: float) -> np.ndarray:
    """Reciprocal inversion for Route A."""
    return a_m / (1.0 - h_eq_m / a_m)


def main() -> None:
    radius_earth_m = 6_371_000.0
    gravity_ref_mps2 = 9.80665
    max_equivalent_height_m = 100_000.0

    h_eq_m = np.linspace(0.0, max_equivalent_height_m, 1001)
    r_b_m = route_b_radius(h_eq_m, radius_earth_m)
    r_a_m = route_a_radius(h_eq_m, radius_earth_m)

    z_b_km = (r_b_m - radius_earth_m) / 1000.0
    z_a_km = (r_a_m - radius_earth_m) / 1000.0
    delta_z_m = r_a_m - r_b_m

    out_dir = Path(__file__).resolve().parent
    out_png = out_dir / "routeA_routeB_radius_compare.png"

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(7.2, 8.0), constrained_layout=True, sharex=True
    )

    ax0.plot(h_eq_m / 1000.0, z_b_km, lw=2.3, label="Route B: cubic inversion")
    ax0.plot(h_eq_m / 1000.0, z_a_km, lw=2.3, label="Route A: reciprocal inversion")
    ax0.set_ylabel("Recovered height z (km)")
    ax0.set_title(
        "Radius inversion under a common equivalent integral height\n"
        f"a = {radius_earth_m/1000:.0f} km, G = {gravity_ref_mps2:.5f} m s$^{{-2}}$"
    )
    ax0.grid(True, alpha=0.25)
    ax0.legend(frameon=False)

    ax1.plot(h_eq_m / 1000.0, delta_z_m, color="black", lw=2.0)
    ax1.axhline(0.0, color="0.6", lw=1.0)
    ax1.set_xlabel("Equivalent integral height h = X / G (km)")
    ax1.set_ylabel("Route A - Route B (m)")
    ax1.grid(True, alpha=0.25)

    h_100 = max_equivalent_height_m
    z_b_100 = route_b_radius(np.array([h_100]), radius_earth_m)[0] - radius_earth_m
    z_a_100 = route_a_radius(np.array([h_100]), radius_earth_m)[0] - radius_earth_m
    diff_100 = z_a_100 - z_b_100

    ax1.annotate(
        f"At h = 100 km:\nΔz = {diff_100:.0f} m",
        xy=(h_100 / 1000.0, diff_100),
        xytext=(70.0, diff_100 * 0.55),
        arrowprops={"arrowstyle": "->", "lw": 1.0},
        fontsize=10,
    )

    fig.savefig(out_png, dpi=220)

    print(f"Saved figure to: {out_png}")
    print(f"Route B recovered height at h=100 km: {z_b_100/1000.0:.6f} km")
    print(f"Route A recovered height at h=100 km: {z_a_100/1000.0:.6f} km")
    print(f"Difference (A - B) at h=100 km: {diff_100:.3f} m")


if __name__ == "__main__":
    main()
