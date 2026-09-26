"""
Validation checks for the Silicon-Photonic WDM Link model.

This script reproduces the numerical validation tests reported in the
project README:

1. Monte Carlo convergence
2. Random-seed robustness
3. Wavelength-grid convergence
4. Analytical vs. numerical radius sensitivity
"""

import numpy as np


# ============================================================
# 1. BASELINE MODEL PARAMETERS
# ============================================================

lambda0 = 1.55          # um
neff0 = 2.151
ng0 = 4.2
R = 10.0               # um

kappa = 0.2
a = 0.98
t = np.sqrt(1 - kappa**2)

sigma_radius_nm = 5.0
sigma_neff = 0.001
max_allowed_shift_nm = 0.5

dneff_dlambda = (neff0 - ng0) / lambda0


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def calculate_beta(wavelength_grid):
    """Calculate propagation constant across a wavelength grid."""

    neff_grid = (
        neff0
        + dneff_dlambda * (wavelength_grid - lambda0)
    )

    beta_grid = 2 * np.pi * neff_grid / wavelength_grid

    return beta_grid


def ring_transmission(R_test, wavelength_grid):
    """Calculate through-port power transmission."""

    beta_grid = calculate_beta(wavelength_grid)

    L_test = 2 * np.pi * R_test
    round_trip_phase = beta_grid * L_test

    transmission = (
        a**2
        + t**2
        - 2 * a * t * np.cos(round_trip_phase)
    ) / (
        1
        + (a * t)**2
        - 2 * a * t * np.cos(round_trip_phase)
    )

    return transmission


def find_resonance_index(wavelength_grid, R_test, m_target):
    """Locate a selected resonance order on a wavelength grid."""

    beta_grid = calculate_beta(wavelength_grid)

    L_test = 2 * np.pi * R_test
    cycles = beta_grid * L_test / (2 * np.pi)

    return np.argmin(np.abs(cycles - m_target))


def calculate_sensitivities():
    """
    Reproduce the numerical radius and effective-index sensitivities
    used in the tolerance model.
    """

    wavelengths = np.linspace(1.50, 1.60, 50000)
    beta = calculate_beta(wavelengths)

    L = 2 * np.pi * R

    m_target = round(
        beta[np.argmin(np.abs(wavelengths - lambda0))]
        * L
        / (2 * np.pi)
    )

    # Radius sensitivity
    radius_errors_nm = np.linspace(-20, 20, 21)

    resonance_shifts_radius_nm = []

    nominal_index = find_resonance_index(
        wavelengths, R, m_target
    )

    nominal_resonance_nm = (
        wavelengths[nominal_index] * 1000
    )

    for radius_error_nm in radius_errors_nm:

        R_test = R + radius_error_nm / 1000

        index = find_resonance_index(
            wavelengths, R_test, m_target
        )

        resonance_nm = wavelengths[index] * 1000

        resonance_shifts_radius_nm.append(
            resonance_nm - nominal_resonance_nm
        )

    radius_fit = np.polyfit(
        radius_errors_nm,
        resonance_shifts_radius_nm,
        1
    )

    radius_sensitivity = radius_fit[0]

    # Effective-index sensitivity
    neff_errors = np.linspace(-0.01, 0.01, 21)

    resonance_shifts_neff_nm = []

    for delta_neff in neff_errors:

        neff_shifted = (
            neff0
            + delta_neff
            + dneff_dlambda * (wavelengths - lambda0)
        )

        beta_shifted = (
            2 * np.pi * neff_shifted / wavelengths
        )

        cycles_shifted = (
            beta_shifted * L / (2 * np.pi)
        )

        index = np.argmin(
            np.abs(cycles_shifted - m_target)
        )

        resonance_nm = wavelengths[index] * 1000

        resonance_shifts_neff_nm.append(
            resonance_nm - nominal_resonance_nm
        )

    neff_fit = np.polyfit(
        neff_errors,
        resonance_shifts_neff_nm,
        1
    )

    neff_sensitivity = neff_fit[0]

    return radius_sensitivity, neff_sensitivity


# ============================================================
# 3. MONTE CARLO FUNCTION
# ============================================================

def run_monte_carlo(
    num_devices,
    seed,
    radius_sensitivity,
    neff_sensitivity
):
    """Run the fabrication-tolerance Monte Carlo simulation."""

    rng = np.random.default_rng(seed)

    radius_errors_nm = rng.normal(
        0,
        sigma_radius_nm,
        num_devices
    )

    neff_errors = rng.normal(
        0,
        sigma_neff,
        num_devices
    )

    radius_shifts_nm = (
        radius_sensitivity * radius_errors_nm
    )

    neff_shifts_nm = (
        neff_sensitivity * neff_errors
    )

    total_shifts_nm = (
        radius_shifts_nm + neff_shifts_nm
    )

    passing = (
        np.abs(total_shifts_nm)
        <= max_allowed_shift_nm
    )

    yield_percent = (
        np.mean(passing) * 100
    )

    shift_std = np.std(total_shifts_nm)

    return yield_percent, shift_std


# ============================================================
# 4. WAVELENGTH-GRID CONVERGENCE
# ============================================================

def wavelength_grid_test(num_samples):
    """Calculate resonance, FSR, FWHM, and Q for a grid size."""

    wavelengths = np.linspace(
        1.50,
        1.60,
        num_samples
    )

    beta = calculate_beta(wavelengths)

    L = 2 * np.pi * R
    cycles = beta * L / (2 * np.pi)

    # Track the resonance order nearest lambda0
    lambda0_index = np.argmin(
        np.abs(wavelengths - lambda0)
    )

    m_target = round(cycles[lambda0_index])

    resonance_index = np.argmin(
        np.abs(cycles - m_target)
    )

    resonance_wavelength = wavelengths[
        resonance_index
    ]

    transmission = ring_transmission(
        R,
        wavelengths
    )

    # --------------------------------------------------------
    # Numerical FSR
    # Find neighboring resonance orders around m_target
    # --------------------------------------------------------

    resonance_wavelengths = []

    for m in [
        m_target - 2,
        m_target - 1,
        m_target,
        m_target + 1,
        m_target + 2
    ]:

        index = np.argmin(
            np.abs(cycles - m)
        )

        resonance_wavelengths.append(
            wavelengths[index]
        )

    resonance_wavelengths = np.array(
        resonance_wavelengths
    )

    resonance_wavelengths = np.sort(
        resonance_wavelengths
    )

    numerical_fsr_nm = (
        np.mean(
            np.diff(resonance_wavelengths)
        )
        * 1000
    )

    # --------------------------------------------------------
    # FWHM/Q
    # Restrict calculation to ±1 nm around target resonance
    # so neighboring microring resonances cannot interfere.
    # --------------------------------------------------------

    window_mask = (
        np.abs(
            wavelengths - resonance_wavelength
        )
        <= 0.001
    )

    local_wavelengths = wavelengths[
        window_mask
    ]

    local_transmission = transmission[
        window_mask
    ]

    local_min_index = np.argmin(
        local_transmission
    )

    min_transmission = local_transmission[
        local_min_index
    ]

    max_transmission = np.max(
        local_transmission
    )

    half_level = (
        min_transmission
        + (max_transmission - min_transmission) / 2
    )

    left_side = local_transmission[
        :local_min_index
    ]

    right_side = local_transmission[
        local_min_index + 1:
    ]

    left_index = np.argmin(
        np.abs(left_side - half_level)
    )

    right_index_local = np.argmin(
        np.abs(right_side - half_level)
    )

    right_index = (
        local_min_index
        + 1
        + right_index_local
    )

    FWHM = (
        local_wavelengths[right_index]
        - local_wavelengths[left_index]
    )

    fwhm_nm = FWHM * 1000

    Q = resonance_wavelength / FWHM

    return (
        numerical_fsr_nm,
        resonance_wavelength * 1000,
        fwhm_nm,
        Q
    )


# ============================================================
# 5. RUN VALIDATION
# ============================================================

radius_sensitivity, neff_sensitivity = (
    calculate_sensitivities()
)

print("\n========================================")
print("MODEL VALIDATION")
print("========================================")

print(
    f"\nNumerical radius sensitivity: "
    f"{radius_sensitivity:.4f} nm/nm"
)

print(
    f"Numerical neff sensitivity: "
    f"{neff_sensitivity:.4f} nm/index-unit"
)


# ------------------------------------------------------------
# Monte Carlo convergence
# ------------------------------------------------------------

print("\n--- Monte Carlo Convergence ---")

for num_devices in [1000, 10000, 100000]:

    yield_percent, shift_std = run_monte_carlo(
        num_devices,
        42,
        radius_sensitivity,
        neff_sensitivity
    )

    print(
        f"N = {num_devices:6d} | "
        f"Yield = {yield_percent:5.2f}% | "
        f"Shift std = {shift_std:.4f} nm"
    )


# ------------------------------------------------------------
# Random-seed robustness
# ------------------------------------------------------------

print("\n--- Random-Seed Robustness ---")

for seed in [1, 10, 100, 1000, 2026]:

    yield_percent, shift_std = run_monte_carlo(
        1000,
        seed,
        radius_sensitivity,
        neff_sensitivity
    )

    print(
        f"Seed = {seed:4d} | "
        f"Yield = {yield_percent:5.2f}% | "
        f"Shift std = {shift_std:.4f} nm"
    )


# ------------------------------------------------------------
# Wavelength-grid convergence
# ------------------------------------------------------------

print("\n--- Wavelength-Grid Convergence ---")

for num_samples in [10000, 50000, 100000]:

    fsr, resonance, fwhm, Q = (
        wavelength_grid_test(num_samples)
    )

    print(
        f"N = {num_samples:6d} | "
        f"FSR = {fsr:.3f} nm | "
        f"Resonance = {resonance:.3f} nm | "
        f"FWHM = {fwhm:.3f} nm | "
        f"Q = {Q:.0f}"
    )


# ------------------------------------------------------------
# Analytical radius-sensitivity check
# ------------------------------------------------------------

analytical_radius_sensitivity = (
    lambda0 * 1000 * neff0
    / ((R * 1000) * ng0)
)

relative_error_percent = (
    abs(
        radius_sensitivity
        - analytical_radius_sensitivity
    )
    / analytical_radius_sensitivity
    * 100
)

print("\n--- Analytical Radius-Sensitivity Check ---")

print(
    f"Numerical:  "
    f"{radius_sensitivity:.5f} nm/nm"
)

print(
    f"Analytical: "
    f"{analytical_radius_sensitivity:.5f} nm/nm"
)

print(
    f"Difference: "
    f"{relative_error_percent:.2f}%"
)

print("\nValidation complete.")