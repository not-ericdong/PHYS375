import numpy as np

import matplotlib.pyplot as plt

##2a,b)
def R_wd(M):
    '''
    Radius of white dwarf equation. Consumes mass of star in units of solar
    radii and returns the log of the radius of the white dwarf in km.
    '''
    R_sun = 696342
    return np.log10(0.01*R_sun*(M/0.7)**(-1.0/3.0))

def R_ns(M):
    '''
    Radius of neutron star equation. Consumes mass of star in units of solar
    radii and returns the log of the radius of the neutron star in km.
    '''
    return np.log10(11*(M/1.4)**(-1.0/3.0))

def R_bh(M):
    '''
    Radius of black hole equation. Consumes mass of BH in units of solar
    radii and returns the log of the radius of the BH in km.
    '''
    return np.log10(3*(M))

def R_conrho(M):
    '''
    Radius of Mass-Radius equation. Consumes mass of star in units of solar
    radii and returns the log of the radius of the star in km that has the
    density of the core of the sun.
    '''
    R_sun = 696342
    return np.log10((3*M/(4*np.pi*(0.599)))**(1.0/3.0)*0.02*R_sun)

M_wd = np.linspace(0.01, 1.4, 1000)
M_ns = np.linspace(1.4, 3, 1000)
M_bh = np.linspace(3, 5, 1000)
M_conrho = np.linspace(0.01, 5, 1000)

R_wd = R_wd(M_wd)
R_ns = R_ns(M_ns)
R_bh = R_bh(M_bh)
R_conrho = R_conrho(M_conrho)

M_wd = np.log10(M_wd)
M_ns = np.log10(M_ns)
M_bh = np.log10(M_bh)
M_conrho = np.log10(M_conrho)

R0, = plt.plot(M_wd,R_wd,label='White Dwarf')
R1, = plt.plot(M_ns,R_ns,label='Neutron Star')
R2, = plt.plot(M_bh,R_bh,label='Black Hole')
R3, = plt.plot(M_conrho, R_conrho, label='M-R Relation')

axes = plt.gca()
R0.set_dashes([10,5])
R1.set_dashes([10,5])

plt.title("Radius and Mass Relation of Stars")

#plt.xlabel('log($\\frac{M}{M_{\\odot}}$)', fontsize=10)
#plt.ylabel('log($\\frac{R}{km}$)', fontsize=10)
lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2, numpoints = 1)
#plt.yticks(fontsize=8)
#plt.xticks(fontsize=8)
plt.savefig('a5q2ba.png', format='png', bbox_extra_artists=(lgd,),bbox_inches='tight')
#plt.show()
plt.close()

##2d)
def w_wd(M):
    '''
    Angular momentum of white dwarf in rev/min with given M in units of
    solar mass
    '''
    return 9.261e-6*M**(-1.0/3.0)

def w_ns(M):
    '''
    Angular momentum of neutron star in rev/min with given M in units of
    solar mass
    '''
    R_sun = 696342
    return 4.8214e-12*M**(-1.0/3.0)*R_sun**2

w_bh = 2.704e-11*696342**2

M_wd = np.linspace(0.01, 1.4, 1000)
M_ns = np.linspace(1.4, 3, 1000)
M_bh = 3

w_wd = w_wd(M_wd)
w_ns = w_ns(M_ns)

w_wd = np.log10(w_wd)
w_ns = np.log10(w_ns)
w_bh = np.log10(w_bh)
M_wd = np.log10(M_wd)
M_ns = np.log10(M_ns)
M_bh = np.log10(M_bh)

R0, = plt.plot(M_wd,w_wd,label='White Dwarf')
R1, = plt.plot(M_ns,w_ns,label='Neutron Star')
plt.scatter(M_bh,w_bh,label='Black Hole')
plt.grid(True)
axes = plt.gca()
R0.set_dashes([4,4])
R1.set_dashes([8,4,2,4])
plt.title('2d) log($\\omega$) vs log($\\frac{M}{M_{\\odot}}$)', fontsize=10)
plt.axis('tight')
plt.xlabel('log($\\frac{M}{M_{\\odot}}$)', fontsize=10)
plt.ylabel('log($\\omega$) (rev/min)', fontsize=10)
lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2, numpoints = 1)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
plt.savefig('a5q2d.png', format='png', bbox_extra_artists=(lgd,),bbox_inches='tight')
#plt.show()
plt.close()
