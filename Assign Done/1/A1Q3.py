import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt

dataraw = ascii.read('W19_assignment1_orbit.dat')
#print(dataraw)

orb_phase = dataraw["col1"]
rad_vel1 = dataraw["col2"]
rad_vel2 = dataraw["col3"]
app_mag = dataraw["col4"]

# part a
time = 50*orb_phase
"""
plt.plot(time, rad_vel1, "b^")
plt.plot(time, rad_vel2, "r^")
plt.xlabel("Time [Days]")
plt.ylabel("Radial Velocity [km/s]")
plt.show()
"""
# part b
vel_sum = max(rad_vel1)*1000 + max(rad_vel2)*1000
#print(vel_sum)
constants = (4.32*10**6.0)/(2*np.pi*6.674*10**(-11))
#print(constants)
masses = constants*vel_sum**3.0

print(masses)
# part c
Lo = (10**(0.4*(min(app_mag))))
#print(Lo)
lume = []
for x in app_mag:
    lume.append(np.log10((10**(0.4*(float(x))))/Lo))
#print(lume)

plt.plot(time, lume, 'g--')
plt.plot(time, lume, 'g^')
plt.xlabel("Time")
plt.ylabel("Luminosity")
#plt.ylim(12.2)
plt.show()

# part d

L1 = max(lume)
print(L1)
print(lume[0])
print(min(lume))

ansr = (L1/lume[0])**0.25
print(ansr)

# part e
count1 = 0
for x in lume:
    count1 += 1
    if x == 0:
        print(time[count1])
        print(rad_vel1[0])
        print(rad_vel2[0])
        print(count1)
        break

count2 = 0
lenx = 0
for x in lume:
    count2 += 1
    if x == max(lume):
        print(time[count2])
        print(rad_vel2[count2])
        print(rad_vel1[count2])
        break





