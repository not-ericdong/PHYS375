import numpy as np
import matplotlib.pyplot as plt

# generating data
temp1 = 3000
temp2 = 5500
temp3 = 30000
h = 6.626*10**(-34)
c = 2.998*10**(8)
c2 = c**2
kb = 1.38*10**(-23)
data_wavelength = []
for i in range(400, 800):
    data_wavelength.append(i)
data_wavelength1 = []
for i in range(50, 3000):
    data_wavelength1.append(i)

# planck function formula in a
b_planck1 = []
b_planck2 = []
b_planck3 = []
def planck(wave, T):
    for i in wave:
        spec = (i*10**-9*((2 * h * c2) / (i*10**-9) ** 5) * \
               (1.0 / (np.exp((h * c) / (i*10**-9 * kb * T)) - 1)))
        if T == 3000:
            b_planck1.append(spec)
        if T == 5500:
            b_planck2.append(spec)
        if T == 30000:
            b_planck3.append(spec)

B_planck2 = np.array(planck(data_wavelength1, temp2))
B_planck3 = np.array(planck(data_wavelength1, temp3))
B_planck1 = np.array(planck(data_wavelength1, temp1))

#plt.plot(data_wavelength1, b_planck1, 'b^')
plt.plot(data_wavelength1, b_planck2, 'r^')
#plt.plot(data_wavelength1, b_planck3, 'y^')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Spectral Radiance (W/m^3*nm)")
#plt.xscale("log")
plt.show()
