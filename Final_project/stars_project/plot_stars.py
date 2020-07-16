import pandas as pd
from math import pi
import numpy as np
import matplotlib.pyplot as plt
import os
sb = 5.67e-8

file_loc = input("Enter file name for Low Mass (i.e. File.csv): ")
data = pd.read_csv(file_loc)

idx = data["idx"][0]
print()
r = data["r"].values
r_s = r[idx]
T = data["T"].values
T_s = T[idx]
T_c = T[0]
rho = data["rho"].values
rho_c = rho[0]
L = data["L"].values
L_s = L[idx]
M = data["M"].values
M_s = M[idx]
print("mass of star:",M_s)
print("temp of star:",T_s)
print("L of star:",L_s)
print("R of star:",r_s)
P = data["P"].values
P_c = P[0]
tau = data["tau"].values
kappa = data["kappa_l"].values
dL_dr = data["dL/dr"]


dlogP=np.diff(np.log10(P))
dlogT=np.diff(np.log10(T))
dlogP_dlogT = np.array(dlogP)/np.array(dlogT)

x = np.array(r[:-1])/r_s

y = np.array(dlogP)/np.array(dlogT)
z = np.diff(y)/np.diff(x)


index = z[int(0.2 * len(z)) : int(0.9* len(z))].argmin()
y_new = y[int(0.2 * len(y)) : int(0.9* len(y))][index]
x_new = x[int(0.2 * len(x)) : int(0.9* len(x))][index]

plt.plot(np.array(r)/(r_s),np.array(T)/T_c)
plt.title("$T/T_s$ vs $R/R_s$")
plt.ylabel("$T/T_s$")
plt.xlabel("$R/R_s$")
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.ylim(0,1.1)
plt.show()

plt.plot(np.array(r)/(r_s),np.array(P)/P_c)
plt.title("$P/P_c$ vs $R/R_s$")
plt.ylabel("$P/P_c$")
plt.xlabel("$R/R_s$")
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.ylim(0,1.1)
plt.show()

plt.plot(np.array(r)/(r_s),np.array(rho)/rho_c)
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.ylim(0,1.1)
plt.title("$\\rho/\\rho_c$ vs $R/R_s$ ")
plt.ylabel("$\\rho/\\rho_c$")
plt.xlabel("$R/R_s$")
plt.show()

plt.plot(np.array(r)/(r_s),np.array(L)/L_s)
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.title("$L/L_s$ vs $R/R_s$")
plt.xlabel("$R/R_s$")
plt.ylabel("$L/L_s$")
plt.show()

plt.plot(np.array(r)/(r_s),np.array(M)/M_s)
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.title("$M/M_s$ vs $R/R_s$")
plt.xlabel("$R/R_s$")
plt.ylabel("$M/M_s$")
plt.show()

plt.plot(np.array(r)/(r_s),np.log10(np.array(tau)))
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.title("$log_{10}\\tau$ vs $R/R_s$")
plt.ylabel("$log10(\\tau)$")
plt.xlabel("$R/R_s$")
plt.show()

plt.plot(np.array(r)/(r_s),np.log10(np.array(kappa)))
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.title("$log_{10}\\kappa$ vs $R/R_s$")
plt.ylabel("$log10(\\kappa)$")
plt.xlabel("$R/R_s$")
plt.show()

###dlopP is 1 value, not an array
plt.plot(np.array(r[:-1])/r_s,dlogP_dlogT)
plt.title("$dlogP/dlogT$ vs $R/R_s$")
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.ylim(0,10)
plt.xlabel("$R/R_s$")
plt.ylabel("$dlog_{10}P/dlog_{10}T$")
plt.show()

dL_dr
plt.plot(np.array(r)/r_s,dL_dr*r_s/L_s)
plt.title("$dL/dr$ vs $R/R_s$")
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.xlabel("$R/R_s$")
plt.ylabel("$dL/dr$")
plt.show()

file_loc = input("Enter file name for High Mass (i.e. File.csv): ")
data = pd.read_csv(file_loc)


idx = data["idx"][0]
print()
r = data["r"].values
r_s = r[idx]
T = data["T"].values
T_s = T[idx]
T_c = T[0]
rho = data["rho"].values
rho_c = rho[0]
L = data["L"].values
L_s = L[idx]
M = data["M"].values
M_s = M[idx]
print("mass of star:",M_s)
print("temp of star:",T_s)
print("L of star:",L_s)
print("R of star:",r_s)
P = data["P"].values
P_c = P[0]
tau = data["tau"].values
kappa = data["kappa_l"].values
dL_dr = data["dL/dr"]


dlogP=np.diff(np.log10(P))
dlogT=np.diff(np.log10(T))
dlogP_dlogT = np.array(dlogP)/np.array(dlogT)

x = np.array(r[:-1])/r_s

y = np.array(dlogP)/np.array(dlogT)
z = np.diff(y)/np.diff(x)


index = z[int(0.6* len(z)) : int(0.99* len(z))].argmin()
y_new = y[int(0.6 * len(y)) : int(0.99* len(y))][index]
x_new = x[int(0.6 * len(x)) : int(0.99* len(x))][index]
index1 = z[int(0.02 * len(z)) : int(0.3* len(z))].argmax()
y_new1= y[int(0.02 * len(y)) : int(0.3* len(y))][index1]
x_new1 = x[int(0.02 * len(x)) : int(0.3* len(x))][index1]

plt.plot(np.array(r)/(r_s),np.array(T)/T_c)
plt.title("$T/T_c$ vs $R/R_s$")
plt.ylabel("$T/T_c$")
plt.xlabel("$R/R_s$")
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.ylim(0,1.1)
plt.show()

plt.plot(np.array(r)/(r_s),np.array(P)/P_c)
plt.title("$P/P_c$ vs $R/R_s$")
plt.ylabel("$P/P_c$")
plt.xlabel("$R/R_s$")
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.ylim(0,1.1)
plt.show()

plt.plot(np.array(r)/(r_s),np.array(rho)/rho_c)
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.ylim(0,1.1)
plt.title("$\\rho/\\rho_c$ vs $R/R_s$ ")
plt.ylabel("$\\rho/\\rho_c$")
plt.xlabel("$R/R_s$")
plt.show()

plt.plot(np.array(r)/(r_s),np.array(L)/L_s)
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.title("$L/L_0$ vs $R/R_s$")
plt.xlabel("$R/R_s$")
plt.ylabel("$L/L_0$")
plt.show()

plt.plot(np.array(r)/(r_s),np.array(M)/M_s)
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.title("$M/M_0$ vs $R/R_s$")
plt.xlabel("$R/R_s$")
plt.ylabel("$M/R_s$")
plt.show()

plt.plot(np.array(r)/(r_s),np.log10(np.array(tau)))
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.title("$log_{10}\\tau$ vs $R/R_s$")
plt.ylabel("$log10(\\tau)$")
plt.xlabel("$R/R_s$")
plt.show()

plt.plot(np.array(r)/(r_s),np.log10(np.array(kappa)))
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.title("$log_{10}\\kappa$ vs $R/R_s$")
plt.ylabel("$log10(\\kappa)$")
plt.xlabel("$R/R_s$")
plt.show()

###dlopP is 1 value, not an array
plt.plot(np.array(r[:-1])/r_s,dlogP_dlogT)
plt.title("$dlogP/dlogT$ vs $R/R_s$")
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.ylim(0,10)
plt.xlabel("$R/R_s$")
plt.ylabel("$dlog_{10}P/dlog_{10}T$")
plt.show()

dL_dr
plt.plot(np.array(r)/r_s,dL_dr*r_s/L_s)
plt.title("$dL/dr$ vs $R/R_s$")
plt.xlim(0,1)
plt.axvspan(x_new,1, color = '0.75')
plt.axvspan(0,x_new1, color = '0.75')
plt.xlabel("$R/R_s$")
plt.ylabel("$dL/dr$")
plt.show()