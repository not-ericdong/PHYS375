import numpy as np
import scipy as sp
from scipy import constants as con
import matplotlib.pyplot as plt

def rho(n):
    x = np.linspace(0,5,10000)
    theta = 1
    dtheta = 0
    ddtheta = -1
    dx = x[1]-x[0]
    theta_list = [1]
    M_list = [0]
    M = 0
    Ps_list = [0]
    x2 = np.linspace(0,5,10000)
    #MJ = 0
    MJ_list = [0]
    x3 = np.linspace(0,4000000000000000,10000)
    for i in range(len(x)-1):
        i += 1
        ddtheta = - (4*np.pi*theta**2 + ((2.0/x[i]) - (1.0/theta) * dtheta)*dtheta)
        dtheta += ddtheta * dx
        theta += dtheta * dx
        theta_list.append(theta)
        dM = 4*np.pi*x[i]**2*theta
        M += dM*dx
        M_list.append(M)
        Ps = theta*M**2*(con.k*10/(2.4*con.m_p))**4*(1.0/con.G)**3*(1.0/(1.989*10**30))**2
        x2[i] = x2[i]*(con.G*1.989*10**30*2.4*con.m_p)/(con.k*10*M*1.496*10**11)
        Ps_list.append(Ps)
        MJ = 0.2*1.989*10**30*((theta*(1.989*10**30/M)*(con.k*10/(con.G*1.989*10**30*2.4*con.m_p))**3)/(3.0*10**15))**(-1.0/2.0)
        MJ_list.append(MJ)
        x3[i] = x3[i]/(1.496*10**11)
    return (theta_list,M_list,Ps_list,x2,MJ_list,x3)

rho0=rho(0)

f0 = rho0[0]
x = np.linspace(0,5,10000)
l0, = plt.plot(x,f0)
axes = plt.gca()
plt.title('${\\rho}$ vs r')
plt.xlabel('r')
plt.ylabel('${\\rho}$')
axes.set_ylim([0,1])
plt.savefig('a4q2b1.png', format='png')
plt.close()

M0 = rho0[1]
l0, = plt.plot(x,M0,label='n=0')
axes = plt.gca()
plt.title('M vs r')
plt.xlabel('r')
plt.ylabel('M')
plt.savefig('a4q2b2.png', format='png')
plt.close()

P0 = rho0[2]
partcd=rho0[3]

l0, = plt.plot(partcd,P0,label='n=0')
axes = plt.gca()
axes.set_xlim([0,75000])
plt.title('Ps vs r')
plt.xlabel('r (AU)')
plt.ylabel('Ps (Pa)')
plt.savefig('a4q2c.png', format='png',)
plt.close()

P0 = rho0[4]
partcd=rho0[3]

l0, = plt.plot(partcd,P0,label='n=0')
axes = plt.gca()
axes.set_xlim([0,50000])
plt.title('Mj vs r')
plt.xlabel('r (AU)')
plt.ylabel('Mj (kg)')
plt.savefig('a4q2d.png', format='png',)
plt.close()