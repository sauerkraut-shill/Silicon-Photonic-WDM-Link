import numpy as np
import matplotlib.pyplot as plt


lambda0 = 1.55
neff0 = 2.151
ng0 = 4.2
wavelengths = np.linspace(1.50, 1.60, 500)

dneff_dlambda = (neff0 - ng0) / lambda0

beta0 = 2 * np.pi * neff0 / lambda0

print("Effective index:", neff0)
print("Propagation constant:", beta0, "1/um")
print("dn_eff/dlambda:", dneff_dlambda, "1/um")

L = 10.0

phase = beta0 * L

print("Accumulated phase:", phase, "radians")
print("Number of optical cycles:", phase / (2 * np.pi))



neff = neff0 + dneff_dlambda * (wavelengths - lambda0)

beta = 2 * np.pi * neff / wavelengths

phase = beta * L

# Phase plot
plt.figure()
plt.plot(wavelengths, phase)
plt.xlabel("Wavelength (um)")
plt.ylabel("Phase (radians)")
plt.title("Accumulated Phase vs Wavelength")
plt.grid()


# Propagation constant plot
plt.figure()
plt.plot(wavelengths, beta)
plt.xlabel("Wavelength (um)")
plt.ylabel("Propagation constant beta (1/um)")
plt.title("Waveguide Propagation Constant vs Wavelength")
plt.grid()


# Effective index plot
plt.figure()
plt.plot(wavelengths, neff)
plt.xlabel("Wavelength (um)")
plt.ylabel("Effective index")
plt.title("Effective Index vs Wavelength")
plt.grid()


# Display all figures
plt.show()