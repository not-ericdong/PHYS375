import numpy as np
import matplotlib.pyplot as plt

# generating data
temp = 5500
h = 6.626*10**(-34)
c = 2.998*10**(8)
c2 = c**2
kb = 1.38*10**(-23)
data_wavelength = []
for i in range(100, 2000):
    data_wavelength.append(i)

# planck function formula in a
b_planck = []
def planck(wave):
    for i in wave:
        spec = ((2 * h * c2) / (i*10**-9) ** 5) * \
               (1.0 / (np.exp((h * c) / (i*10**-9 * kb * temp)) - 1))
        if spec == 20563804981614.773:
            print("this is the max", i)
        b_planck.append(spec)

B_planck = np.array(planck(data_wavelength))
print(b_planck)
print(max(b_planck))

plt.plot(data_wavelength, b_planck, 'b^')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Spectral Radiance (W/m^3)")
plt.show()

