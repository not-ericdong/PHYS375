import numpy as np
import matplotlib.pyplot as plt

# generating data

p = 10**-6
k = 3
Te = 10000
data_distance = []
for i in range(0, 1000*10**3):
    data_distance.append(i)

# planck function formula, Bv
yaxisdata = []
def function(a):
    for i in a:
        edding = ((3.0/4.0)*Te**4*(k*p*i+(2.0/3.0)))**(1.0/4.0)
        if edding >= 10000 and edding <= 10005:
            print("hi", i)
        yaxisdata.append(edding)

arrayray = np.array(function(data_distance))
#print(yaxisdata)

plt.plot(data_distance, yaxisdata, 'b^')
plt.xlabel("Distance (m)")
plt.ylabel("Temperature (K)")
#plt.xlim(10**8, 3.1*10**15)
#plt.gca().invert_yaxis()
plt.show()
