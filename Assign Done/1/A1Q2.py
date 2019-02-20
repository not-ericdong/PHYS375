import matplotlib.pyplot as plt
import numpy as np

# open files #
data = open("hipparcos.txt","r")
graphlabel = open("hipparcos_cols.txt", "r")

# sort data into columns #
lines = data.readlines()
data_column1 = []   # parallax in (milliarcseconds)
data_column2 = []   # V magnitude
data_column3 = []   # B magnitude
data_column4 = []   # I magnitude
for x in lines:
    data_column1.append(x.split()[0])
    data_column2.append(x.split()[1])
    data_column3.append(x.split()[2])
    data_column4.append(x.split()[3])

graphlabel_list = []
graph_lines = graphlabel.readlines()
for x in graph_lines:
    graphlabel_list.append(x)
# print(graphlabel_list)

# ---Calculations--- #
# Get B-V #
data_BV = []
for star in lines:
    #B mags - V mags
    data_BV.append(float(star.split()[2]) - float(star.split()[1]))
# print(data_BV)

# Get absolute V magnitude #
Mv = []
for star in lines:
    #M = m-5log10(distance in parsecs)-1 and 1/parallax is distance in parsecs - 1
    Mv.append(float(star.split()[1])-5*(np.log10(1.0/(float(star.split()[0])*0.001)))-1)
# print(Mv)

# Get Temperature #
temp = []
for star in data_BV:
    #equation to find temp
    temp.append(9000/(star+0.93))
# print(temp)

# Get Lv/Lsol #
lume = []
for star in Mv:
    #10^-0.4(Mv-Msol)
    lume.append(10 ** (-0.4*(float(star) - 4.83)))
# print(lume)

# Stefan-Boltzmann_law #
sun_rad = 695508000
lumesun = 3.828*10**26.0
SBconst = 5.67*10.0**(-8.0)
const = 4*np.pi*SBconst
def stefb(R,T):
    return const*R**2.0*T**4.0

# close files #
data.close()
graphlabel.close()

# graphing time #
#2a
x1 = np.array(data_BV)
y1 = np.array(Mv)
plt.plot(x1, y1, "bs")
plt.gca().invert_yaxis()
plt.xlabel("(B-V) Colour")
plt.ylabel("Absolute V Magnitude")
plt.show()
#2b&c
x2 = np.log10(np.array(temp))
y2 = np.log10(np.array(lume))
line1 = []
line2 = []
line3 = []
liney = []

for x in range(len(temp)):
    one = np.log10(stefb(sun_rad, temp[x])/lumesun)
    two = np.log10(stefb(0.2*sun_rad, temp[x])/lumesun)
    three = np.log10(stefb(5*sun_rad, temp[x])/lumesun)
    line1.append(one)
    line2.append(two)
    line3.append(three)

x3 = np.log10(temp)

# lines in part c
plt.plot(x3, line1, dashes=[5,1,10,1], color="blue")
plt.plot(x3, line2, dashes=[2,5,2,10], color="red")
plt.plot(x3, line3, dashes=[10,10,1,1], color="orange")
# plot for b
plt.plot(x2, y2, "gs")
plt.xlabel("Temperature [K]")
plt.ylabel("Luminosity")
#plt.yscale('log')
#plt.xscale('log')
plt.show()


#fig1 = plt.figure(figsize=[8,6])
#plt.savefig(fig1,dpi=300,bbox_inches="tight")
