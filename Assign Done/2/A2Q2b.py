import numpy as np
import matplotlib.pyplot as plt

# generating data
temp = 5500
h = 6.626*10**(-34)
c = 2.998*10**(8)
c2 = c**2
kb = 1.38*10**(-23)
data_wavelength = []
for i in range(15*10**12, 3000*10**12, 10**10):
    data_wavelength.append(i)

# planck function formula, Bv
b_planck = []
def planck(wave):
    for i in wave:
        spec = ((2*h*i**3)/(c2))*\
               (1.0/(np.exp((h*i)/(kb*temp))-1))
        b_planck.append(spec)

B_planck = np.array(planck(data_wavelength))
print(b_planck)

plt.plot(data_wavelength, b_planck, 'b^')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Spectral Radiance (W/m^2*Hz)")
#plt.xlim(10**8, 3.1*10**15)
#plt.gca().invert_yaxis()
plt.show()

