"""
Silicon Photonic WDM Microring Filter Model

Models the through-port response of dispersive silicon microring resonators
and uses radius tuning to construct a four-channel WDM filter bank.

The model includes:
- Wavelength-dependent effective index
- Microring resonance and free spectral range (FSR)
- Resonance linewidth, Q factor, and extinction ratio
- Resonance sensitivity to ring radius
- Four-channel WDM filter-bank design
- On-channel and off-channel transmission analysis

Wavelengths and geometric dimensions are calculated internally in micrometers
unless otherwise specified. Results are primarily reported in nanometers.
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1. Physical Parameters and Wavelength Grid
# ============================================================

# Nominal optical/device parameters
lambda0 = 1.55
neff0 = 2.151
ng0 = 4.2
R = 10.0
kappa = 0.2 #Coupling coefficient
a = 0.98 # a is the field-amplitude transmission factor after one round trip.

# WDM design parameter
channel_spacing_nm = 2.0

# Fabrication-analysis assumptions
num_devices = 1000
sigma_radius_nm = 5.0
sigma_neff = 0.001
max_allowed_shift_nm = 0.5

# Thermal-model assumptions
reference_temperature_C = 25.0
thermo_optic_coefficient = 1.86e-4


L = 2 * np.pi * R

wavelengths = np.linspace(1.50, 1.60, 50000)

# First-order effective-index dispersion model.
# Using ng = neff - λ(dneff/dλ), estimate the slope at lambda0
# and use it to approximate neff across the wavelength grid.

dneff_dlambda = (neff0 - ng0) / lambda0

neff = neff0 + dneff_dlambda * (wavelengths - lambda0)

beta = 2 * np.pi * neff / wavelengths


t = np.sqrt(1 - kappa**2)  # kappa and t are field coupling coefficients satisfying t² + kappa² = 1.


def ring_transmission(R, wavelength_grid=None):
    """
    Calculate the through-port power transmission of an all-pass microring.

    Parameters
    ----------
    R : float
        Ring radius in micrometers.
    wavelength_grid : ndarray, optional
        Wavelength grid in micrometers. If not provided, use the
        baseline global wavelength grid.

    Returns
    -------
    transmission : ndarray
        Through-port power transmission across the wavelength grid.
    """

    if wavelength_grid is None:
        wavelength_grid = wavelengths

    neff_grid = neff0 + dneff_dlambda * (wavelength_grid - lambda0)
    beta_grid = 2 * np.pi * neff_grid / wavelength_grid

    L = 2 * np.pi * R
    round_trip_phase = beta_grid * L

    transmission = (
        a**2 + t**2 - 2 * a * t * np.cos(round_trip_phase)
    ) / (
        1 + (a * t)**2 - 2 * a * t * np.cos(round_trip_phase)
    )

    return transmission

def find_resonance_wavelength(R_test, m_target):
    """
    Find the wavelength corresponding to a selected resonance order.

    Parameters
    ----------
    R_test : float
        Ring radius in micrometers.
    m_target : int
        Resonance order to track.

    Returns
    -------
    resonance_wavelength : float
        Resonance wavelength in micrometers.
    """

    L_test = 2 * np.pi * R_test

    phase_test = beta * L_test

    cycles_test = phase_test / (2 * np.pi)

    index = np.argmin(np.abs(cycles_test - m_target))

    resonance_wavelength = wavelengths[index]

    return resonance_wavelength

def calculate_resonance_metrics(wavelengths, transmission, resonance_index, R):

    resonance_wavelength = wavelengths[resonance_index]

    fine_wavelengths = np.linspace(
    resonance_wavelength - 0.00001,
    resonance_wavelength + 0.00001,
    10001
    )

    fine_transmission = ring_transmission(R, fine_wavelengths)
    min_transmission = np.min(fine_transmission)
    max_transmission = np.max(transmission)
    half_level = min_transmission + (max_transmission - min_transmission) / 2

    left_side = transmission[:resonance_index]
    left_index = np.argmin(np.abs(left_side - half_level))

    right_side = transmission[resonance_index + 1:]
    right_index_local = np.argmin(np.abs(right_side - half_level))
    right_index = resonance_index + 1 + right_index_local

    FWHM = wavelengths[right_index] - wavelengths[left_index]
    fwhm_nm = FWHM * 1000

    Q = wavelengths[resonance_index] / FWHM

    extinction_ratio_db = 10 * np.log10(
    max_transmission / min_transmission
    )

    return fwhm_nm, Q, extinction_ratio_db, half_level


def calculate_yield(resonance_shifts_nm, max_allowed_shift_nm):
    passing = np.abs(resonance_shifts_nm) <= max_allowed_shift_nm
    yield_percent = np.mean(passing) * 100
    return yield_percent

print(f"Ring radius: {R:.3f} um")
print(f"Ring circumference: {L:.3f} um")

# ============================================================
# 2. Ring Resonance Analysis
# ============================================================

m0 = neff0 * L / lambda0

print(f"Approximate resonance order: {m0:.2f}")

fsr = lambda0**2 / (ng0 * L)

print(f"Approximate FSR: {fsr * 1000:.3f} nm")

round_trip_phase = beta * L

round_trip_cycles = round_trip_phase / (2 * np.pi)

resonance_orders = np.arange(
    np.ceil(round_trip_cycles.min()),
    np.floor(round_trip_cycles.max()) + 1
)

resonance_wavelengths = []

for m in resonance_orders:
    index = np.argmin(np.abs(round_trip_cycles - m))
    resonance_wavelengths.append(wavelengths[index])

resonance_wavelengths = np.array(resonance_wavelengths)

simulated_fsr = np.abs(np.diff(resonance_wavelengths)) * 1000

fsr_error_percent = (
    abs(np.mean(simulated_fsr) - fsr * 1000)
    / (fsr * 1000)
    * 100
)

print(f"Mean simulated FSR: {np.mean(simulated_fsr):.3f} nm")

print(f"FSR model error: {fsr_error_percent:.3f}%")

transmission = ring_transmission(R)

# ============================================================
# 3.  Resonance Characterization
# ============================================================


# Find the known resonance closest to lambda0
target_resonance = resonance_wavelengths[
    np.argmin(np.abs(resonance_wavelengths - lambda0))
]

# Use a physical wavelength window of +/- 1 nm
window_width = 0.001  # um = 1 nm

mask = (
    (wavelengths >= target_resonance - window_width)
    & (wavelengths <= target_resonance + window_width)
)

local_wavelengths = wavelengths[mask]
local_transmission = transmission[mask]

# Find minimum transmission in this region
local_min_index = np.argmin(local_transmission)

resonance_wavelength = local_wavelengths[local_min_index]
T_min = local_transmission[local_min_index]

print("\nResonance characterization:")
print(f"Resonance wavelength: {resonance_wavelength * 1000:.3f} nm")
print(f"Minimum transmission: {T_min:.6f}")

fwhm_nm, Q, extinction_ratio_db, half_level = calculate_resonance_metrics(
    local_wavelengths,
    local_transmission,
    local_min_index,
    R
)

fwhm = fwhm_nm / 1000

print(f"FWHM: {fwhm_nm:.3f} nm")
print(f"Q factor: {Q:.0f}")
print(f"Extinction ratio: {extinction_ratio_db:.2f} dB")

# ============================================================
# 4. Radius Sensitivity Analysis
# ============================================================


radius_values = np.linspace(9.9, 10.1, 21)

shifted_resonances = []

m_target = round(
    neff0 * L / resonance_wavelength
)

for R_test in radius_values:
    resonance = find_resonance_wavelength(R_test, m_target)
    shifted_resonances.append(resonance * 1000)


shifted_resonances = np.array(shifted_resonances)

# ============================================================
# 5. Four-Channel WDM Filter Bank
# ============================================================

channel_spacing_um = channel_spacing_nm / 1000

# Estimate the radius change needed for the desired channel spacing.
# From the first-order resonance sensitivity:
# Δλ / λ ≈ (neff / ng) * (ΔR / R)
# Therefore:
# ΔR ≈ R * (ng / neff) * (Δλ / λ)

radius_spacing = (
    R
    * (ng0 / neff0)
    * (channel_spacing_um / lambda0)
)

radius_spacing = (
    R
    * (ng0 / neff0)
    * (channel_spacing_um / lambda0)
)

ring_radii = [
    R,
    R + radius_spacing,
    R + 2 * radius_spacing,
    R + 3 * radius_spacing
]

transmissions = []

for radius in ring_radii:
    T = ring_transmission(radius)
    transmissions.append(T)

wdm_mask = (
    (wavelengths >= 1.549)
    & (wavelengths <= 1.559)
)

wdm_wavelengths = wavelengths[wdm_mask]

resonance_wavelengths_wdm = []

for T in transmissions:
    T_local = T[wdm_mask]
    min_index = np.argmin(T_local)
    lambda_res = wdm_wavelengths[min_index]
    resonance_wavelengths_wdm.append(lambda_res * 1000)

print("\nWDM filter bank:")

for i, (radius, lambda_res) in enumerate(
    zip(ring_radii, resonance_wavelengths_wdm)
):
    print(
        f"Ring {i + 1}: "
        f"R = {radius:.4f} um, "
        f"resonance = {lambda_res:.3f} nm"
    )

radius_tuning_valid = np.all(
    np.diff(resonance_wavelengths_wdm) > 0
)

print(f"Radius tuning validation: {radius_tuning_valid}")

channel_spacings = np.diff(resonance_wavelengths_wdm)

print("Channel spacings:")

for i, spacing in enumerate(channel_spacings):
    print(
        f"Channel {i + 1} -> {i + 2}: "
        f"{spacing:.3f} nm"
    )

print(f"Mean channel spacing: {np.mean(channel_spacings):.3f} nm")

spacing_to_linewidth_ratio = (
    np.mean(channel_spacings)
    / (fwhm * 1000)
)

print(
    f"Channel spacing / FWHM: "
    f"{spacing_to_linewidth_ratio:.2f}"
)

wdm_bandwidth = (
    resonance_wavelengths_wdm[-1]
    - resonance_wavelengths_wdm[0]
)

fsr_nm = np.mean(simulated_fsr)

wdm_within_fsr = wdm_bandwidth < fsr_nm

print(f"WDM bandwidth: {wdm_bandwidth:.3f} nm")
print(f"WDM bandwidth within one FSR: {wdm_within_fsr}")

# ============================================================
# 6. Channel Selectivity Analysis
# ============================================================


channel_indices = []

for lambda_channel_nm in resonance_wavelengths_wdm:

    lambda_channel_um = lambda_channel_nm / 1000

    index = np.argmin(
        np.abs(wavelengths - lambda_channel_um)
    )

    channel_indices.append(index)

channel_transmission_matrix = []

for T in transmissions:

    ring_channel_transmissions = []

    for index in channel_indices:

        transmission_value = T[index]

        ring_channel_transmissions.append(transmission_value)

    channel_transmission_matrix.append(ring_channel_transmissions)


channel_transmission_matrix = np.array(channel_transmission_matrix)

on_channel_transmissions = []

for R_ring, lambda_channel_nm in zip(ring_radii, resonance_wavelengths_wdm):

    lambda_channel_um = lambda_channel_nm / 1000

    fine_wavelengths = np.linspace(
        lambda_channel_um - 0.00001,
        lambda_channel_um + 0.00001,
        10001
    )

    fine_transmission = ring_transmission(R_ring, fine_wavelengths)

    on_channel_transmissions.append(
        np.min(fine_transmission)
    )

on_channel_transmissions = np.array(on_channel_transmissions)

off_diagonal_mask = ~np.eye(
    len(transmissions),
    dtype=bool
)

off_channel_transmissions = channel_transmission_matrix[off_diagonal_mask]

worst_off_channel_transmission = np.min(off_channel_transmissions)

worst_off_channel_penalty = 1 - worst_off_channel_transmission

worst_off_channel_penalty_percent = worst_off_channel_penalty * 100

worst_off_channel_transmission_db = 10 * np.log10(worst_off_channel_transmission)

on_channel_transmissions_db = 10 * np.log10(on_channel_transmissions)

worst_on_channel_transmission_db = np.max(on_channel_transmissions_db)

print("\nChannel selectivity:")

for i, transmission_db in enumerate(on_channel_transmissions_db):
    print(
        f"Ring {i + 1} on-channel transmission: "
        f"{transmission_db:.2f} dB"
    )

print(
    f"Worst on-channel transmission: "
    f"{worst_on_channel_transmission_db:.2f} dB"
)

print(
    f"Worst off-channel transmission: "
    f"{worst_off_channel_transmission_db:.4f} dB"
)

print(
    f"Worst off-channel power penalty: "
    f"{worst_off_channel_penalty_percent:.3f}%"
)


# ============================================================
# 7. Radius Fabrication Tolerance
# ============================================================

radius_errors_nm = np.linspace(-20, 20, 21)

radius_errors_um = radius_errors_nm / 1000

fabricated_radii = R + radius_errors_um

fabricated_resonances_nm = []

m_target = round(
    neff0 * L / resonance_wavelength
)

for R_fab in fabricated_radii:
    resonance = find_resonance_wavelength(R_fab, m_target)
    resonance = resonance * 1000
    fabricated_resonances_nm.append(resonance)

fabricated_resonances_nm = np.array(fabricated_resonances_nm)

resonance_shifts_nm = fabricated_resonances_nm - (resonance_wavelength * 1000)

slope, intercept = np.polyfit(
    radius_errors_nm,
    resonance_shifts_nm,
    1
)

print(f"Radius sensitivity: {slope:.4f} nm/nm")
print(f"Fit intercept: {intercept:.4f} nm")

# ============================================================
# 8. Effective Index Fabrication Tolerance
# ============================================================

neff_errors = np.linspace(-0.01, 0.01, 21)

neff_resonances_nm = []

for neff_error in neff_errors:
    neff_fab = neff + neff_error

    beta_fab = (neff_fab * 2 * np.pi) / wavelengths

    phase_fab = beta_fab * L

    cycles_fab = phase_fab / (2 * np.pi)

    index = np.argmin(np.abs(cycles_fab - m_target))

    neff_resonances_nm.append(wavelengths[index] * 1000)

neff_resonances_nm = np.array(neff_resonances_nm)

neff_resonance_shifts_nm = neff_resonances_nm - (resonance_wavelength * 1000)

neff_slope, neff_intercept = np.polyfit(
    neff_errors,
    neff_resonance_shifts_nm,
    1
)

print(f"Effective-index sensitivity: {neff_slope:.4f} nm/index-unit")
print(f"Fit intercept: {neff_intercept:.4f} nm")

neff_fitted_shifts_nm = neff_slope * neff_errors + neff_intercept

# ============================================================
# 9. Monte Carlo Fabrication Analysis
# ============================================================

rng = np.random.default_rng(42)

mc_radius_errors_nm = rng.normal(
    0,
    sigma_radius_nm,
    num_devices
)

print("Radius error mean:", np.mean(mc_radius_errors_nm))
print("Radius error std:", np.std(mc_radius_errors_nm))

mc_radius_shifts_nm = slope * mc_radius_errors_nm

mc_radius_shift_mean = np.mean(mc_radius_shifts_nm)
mc_radius_shift_std = np.std(mc_radius_shifts_nm)

print("Radius-induced resonance shift mean:", mc_radius_shift_mean)
print("Radius-induced resonance shift std:", mc_radius_shift_std)

mc_neff_errors = rng.normal(0, sigma_neff, num_devices)

mc_neff_shifts_nm = neff_slope * mc_neff_errors

mc_neff_shift_mean = np.mean(mc_neff_shifts_nm)
mc_neff_shift_std = np.std(mc_neff_shifts_nm)

print("Index-induced resonance shift mean:", mc_neff_shift_mean)
print("Index-induced resonance shift std:", mc_neff_shift_std)

mc_total_shifts_nm = mc_radius_shifts_nm + mc_neff_shifts_nm

mc_total_shift_mean = np.mean(mc_total_shifts_nm)
mc_total_shift_std = np.std(mc_total_shifts_nm)

print("Total resonance shift mean:", mc_total_shift_mean)
print("Total resonance shift std:", mc_total_shift_std)

# ============================================================
# 10. WDM Fabrication Yield
# ============================================================

passing_devices = np.abs(mc_total_shifts_nm) <= max_allowed_shift_nm
num_passing = np.sum(passing_devices)

yield_percent = calculate_yield(
    mc_total_shifts_nm,
    max_allowed_shift_nm
)

print("Passing devices:", num_passing, "/", num_devices)
print(f"Fabrication yield: {yield_percent:.1f}%")

shift_tolerances_nm = np.linspace(0.1, 1.0, 10)

yield_vs_tolerance = []

for tolerance in shift_tolerances_nm:
    yield_current = calculate_yield(
        mc_total_shifts_nm,
        tolerance
    )
    yield_vs_tolerance.append(yield_current)

radius_variance = mc_radius_shift_std ** 2
neff_variance = mc_neff_shift_std ** 2

print("Radius variance:", radius_variance)
print("Effective-index variance:", neff_variance)

total_variance = radius_variance + neff_variance

radius_contribution_percent = (radius_variance / total_variance) * 100
neff_contribution_percent = (neff_variance / total_variance) * 100

print(f"Radius contribution: {radius_contribution_percent:.1f}%")
print(f"Effective-index contribution: {neff_contribution_percent:.1f}%")

sources = ["Radius Variation", "Effective-Index Variation"]
variance_contributions = [radius_contribution_percent, neff_contribution_percent]

# ============================================================
# 11. Temperature Dependence
# ============================================================

temperatures_C = np.linspace(0, 80, 81)

temperature_changes_C = temperatures_C - reference_temperature_C

thermal_neff_changes = temperature_changes_C * thermo_optic_coefficient

thermal_resonance_shifts_nm = neff_slope * thermal_neff_changes

print("Shift at 0 C:", thermal_resonance_shifts_nm[0])
print("Shift at 25 C:", thermal_resonance_shifts_nm[25])
print("Shift at 80 C:", thermal_resonance_shifts_nm[-1])

thermal_sensitivity_nm_per_C = neff_slope * thermo_optic_coefficient

print(f"Thermal sensitivity: {thermal_sensitivity_nm_per_C:.4f} nm/C")


max_temperature_change_C = max_allowed_shift_nm / thermal_sensitivity_nm_per_C

print(f"Maximum temperature change before exceeding spec: {max_temperature_change_C:.2f} C")


# ============================================================
# 12. Combined Fabrication and Thermal Robustness
# ============================================================

combined_yield_vs_temperature = []

for thermal_shift in thermal_resonance_shifts_nm:
    combined_shifts_nm = thermal_shift + mc_total_shifts_nm
    
    yield_at_temperature = calculate_yield(
        combined_shifts_nm,
        max_allowed_shift_nm
    )
    
    combined_yield_vs_temperature.append(yield_at_temperature)

print("Yield at 0 C:", combined_yield_vs_temperature[0])
print("Yield at 25 C:", combined_yield_vs_temperature[25])
print("Yield at 80 C:", combined_yield_vs_temperature[-1])

combined_yield_vs_temperature = np.array(combined_yield_vs_temperature)

max_yield_index = np.argmax(combined_yield_vs_temperature)

max_yield = combined_yield_vs_temperature[max_yield_index]

max_yield_temperature_C = temperatures_C[max_yield_index]

print("Maximum yield:", max_yield)
print("Temperature at maximum yield:", max_yield_temperature_C)

print(f"Maximum combined yield: {max_yield:.1f}%")
print(f"Temperature at maximum yield: {max_yield_temperature_C:.1f} C")

# ============================================================
# 13. Visualization
# ============================================================


# Plot 1: Ring resonance condition
plt.figure()

plt.plot(
    wavelengths * 1000,
    round_trip_cycles,
    label="Round-trip phase"
)

plt.scatter(
    resonance_wavelengths * 1000,
    resonance_orders,
    label="Resonances"
)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Round-trip phase cycles")
plt.title("Microring Resonance Condition")
plt.grid()
plt.legend()


# Plot 2 LEFT:  Ring Transmission Spectrum

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(wavelengths * 1000, transmission)

ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel("Through-port power transmission")
ax1.set_title("(a) Full Transmission Spectrum")
ax1.grid(True)

# Plot 2 RIGHT: Resonance characterization

ax2.plot(
    local_wavelengths * 1000,
    local_transmission,
    label="Transmission"
)

ax2.axhline(
    half_level,
    linestyle="--",
    label="Half depth level"
)

ax2.axvline(
    resonance_wavelength * 1000,
    linestyle=":",
    label="Resonance"
)

ax2.set_xlim(
    resonance_wavelength * 1000 - 1,
    resonance_wavelength * 1000 + 1
)

ax2.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("Through-port power transmission")
ax2.set_title("(b) Resonance Characterization")
ax2.grid(True)
ax2.legend()

fig.suptitle("Microring Resonator Spectral Characterization")

plt.tight_layout()


# Plot 4: Radius sensitivity
plt.figure()

plt.plot(
    radius_values,
    shifted_resonances,
    marker="o"
)

plt.xlabel("Ring radius (µm)")
plt.ylabel("Resonance wavelength (nm)")
plt.title("Resonance Sensitivity to Ring Radius")
plt.grid()

# Plot 5: Four-channel WDM filter bank

plt.figure()

for i, (T, lambda_res) in enumerate(
    zip(transmissions, resonance_wavelengths_wdm)
):
    plt.plot(
        wavelengths * 1000,
        T,
        label=f"Channel {i + 1}: {lambda_res:.1f} nm"
    )

plt.xlim(1550, 1559)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Through-port power transmission")
plt.title("Four-Channel Microring WDM Filter Bank")
plt.grid()
plt.legend()

#Plot 6 LEFT: Radius Errors Influence on Resonance Shift

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(
    radius_errors_nm,
    resonance_shifts_nm,
    "o",
    label="Simulated Data"
)

ax1.plot(
    radius_errors_nm,
    resonance_shifts_nm,
    label="Linear Fit"
)

ax1.set_xlabel("Radius Error (nm)")
ax1.set_ylabel("Resonance Shift (nm)")
ax1.set_title("(a) Radius Fabrication Error")
ax1.grid(True)
ax1.legend()

#Plot 6 RIGHT: Resonance Sensitivity to Effective-Index Error

ax2.plot(
    neff_errors,
    neff_resonance_shifts_nm,
    "o",
    label="Simulated Data"
)

ax2.plot(
    neff_errors,
    neff_fitted_shifts_nm,
    label="Linear Fit"
)

ax2.set_xlabel("Effective Index Error")
ax2.set_ylabel("Resonance Shift (nm)")
ax2.set_title("(b) Effective-Index Error")
ax2.grid(True)
ax2.legend()

fig.suptitle("Fabrication Sensitivity Analysis")

plt.tight_layout()

#Plot 7 LEFT: Monte Carlo Distribution of Resonance Shift

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.hist(
    mc_total_shifts_nm,
    bins=30
)

ax1.axvline(
    mc_total_shift_mean,
    linestyle="--",
    label="Mean"
)

ax1.axvline(
    mc_total_shift_mean + mc_total_shift_std,
    linestyle="--",
    label="+1σ"
)

ax1.axvline(
    mc_total_shift_mean - mc_total_shift_std,
    linestyle="--",
    label="-1σ"
)

ax1.set_xlabel("Total Resonance Shift (nm)")
ax1.set_ylabel("Number of Devices")
ax1.set_title("(a) Monte Carlo Resonance Shift")
ax1.grid(True)
ax1.legend()

#Plot 7 RIGHT: Resonance Shift Tolerance vs. Fabrication Yield

ax2.plot(shift_tolerances_nm, yield_vs_tolerance)

ax2.set_xlabel("Allowable Resonance Shift (nm)")
ax2.set_ylabel("Fabrication Yield (%)")
ax2.set_title("(b) Fabrication Yield vs. Tolerance")

ax2.scatter(0.5, 66.4, label="Chosen Spec (±0.5 nm)")
ax2.legend()

ax2.grid(True)

fig.suptitle("Fabrication Tolerance Analysis")

plt.tight_layout()

#Plot 8 LEFT: Temperature-Induced Resonance Shift

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(
    temperatures_C,
    thermal_resonance_shifts_nm,
    label="Thermal Resonance Shift"
)

ax1.axhline(max_allowed_shift_nm,  linestyle="--", label="Upper Spec Limit")
ax1.axhline(-max_allowed_shift_nm,  linestyle="--", label="Lower Spec Limit")

ax1.axvline(
    reference_temperature_C,
    linestyle=":",
    label="Reference Temperature"
)

ax1.set_xlabel("Temperature (°C)")
ax1.set_ylabel("Resonance Shift (nm)")
ax1.set_title("(a) Temperature-Induced Resonance Shift")
ax1.grid(True)
ax1.legend()

#Plot 8 RIGHT: Combined Fabrication and Thermal Yield vs. Temperature

ax2.plot(temperatures_C, combined_yield_vs_temperature)

ax2.axvline(
    reference_temperature_C,
    linestyle=":",
    label="Reference Temperature"
)

ax2.set_xlabel("Temperature (°C)")
ax2.set_ylabel("Fabrication Yield (%)")
ax2.set_title("(b) Combined Yield vs. Temperature")
ax2.grid(True)
ax2.legend()

fig.suptitle("Thermal Robustness Analysis")

plt.tight_layout()

# Display all figures at the same time
plt.show()

