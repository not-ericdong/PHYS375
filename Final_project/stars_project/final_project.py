# ---------- Table of Contents ----------
"""
00. Import necessary modules
01. Import scientific constants
02. Metallicities [X, Y, Z]
    02a. Solar Metallicity                                      [X, Y, Z]
	02b. Zero Metallicity Function                              Metal_Zero(X) (Output: [X, Y, 0])
03. Central Boundary Conditions ---
04. Total Pressure, P (Equation 5)                              Total_P(rho, T, X, Y, Z)
05. Mean particle weight (Equation 6)                           mu(X, Y, Z):
06. dP/d(rho) (Equation 7a)                                     dP_drho(rho, T, X, Y, Z)
07. dP/dT (Equation 7b)                                         dP_dT(rho, T, X, Y, Z)
08. Specific Energy Generation Rates (Equations 8a, 8b, 9)
    08a. Proton-Proton Chain (Equation 8a)                      epsilon_pp(rho, T, X)
	08b. CNO (Equation 8b)                                      epsilon_CNO(rho, T, X)
	08c. Total Specific Energy Generation Rate (Equation 9)     epsilon(rho, T, X)
09. Rosseland Mean Opacities (Equations 10, 11, 12)
                                                                Kappa_es(X)
                                                                Kappa_ff(Z, rho, T)
                                                                Kappa_H_minus(Z, rho, T)
10. Radiative Opacity (Equation 14)                             Kappa_Rad(rho, T, X, Z)
11. Differential Equations of Solve (Equation 2)
	11a.  Density Equations (2a)                                drho_dr(M, r, rho, dP_dT, dT_dr, dP_drho)
	11b. Temperature Equations (2b)                             dT_dr(r, T)
	11f. Mass Equations (Equation 2c)                           dM_dr(r, rho)
	11g. Luminosity equation (Equation 2d)                      dL_dr(r, rho, epsilon)
	11h. Optical depth equation (Equation 2e)                   dtau_dr(kappa, rho):
12. Runge-Kutta 4th Order DE Solver                             rk4(Equation, IC_x_0, IC_y_0, Final_x, Resolution)

"""

# ---------- Define metallicities ----------
"""
Metallicities here are defined as [X, Y, Z], where:
    - X = hydrogen fraction
    - Y = Helium fraction
    - Z = Metallicity
Solar values taken from https://en.wikipedia.org/wiki/Metallicity
"""

# ----- Solar Metallicity -----
Metal_Sol = [0.7381, 0.2485, 0.0134]

# ----- Zero Metallicity Function -----
"""
This function is defined for convenience. Input the hydrogen fraction, and it
will return the metallicity [X, Y, 0]
"""


def Metal_Zero(X):
    return [X, 1 - X, 0]


# ----- Equations to be solved -----

# --- Central Boundary Conditions ---
"""
The project description mentions the 'shooting' method (page 4), of selecting
0 for initial values.

NOTE: critical values of rho (rho_c) and temperature (T_c) must be chosen.
"""

M_0 = 0
L_0 = 0

# ---------- Import necessary modules -----------
import matplotlib.pyplot as plt
from numpy import pi, arange, sqrt
import numpy as np
from math import isnan
from scipy.optimize import bisect
from time import time
from numpy import column_stack
import pandas as pd
# ---------- Import scientific constants ----------
from scipy import constants as Co

G = Co.G
hbar = Co.hbar
k_b = Co.k
m_p = Co.m_p
m_e = Co.m_e
c = Co.c
a = 4 * Co.sigma / Co.c
sigma = Co.sigma
gamma = 5 / 3.0
DEBUG_FLAG = 0

# --- Varying Parameters ---
"""
change T_c to generate different stars for your specified metalicity 
if the star is not converging increase the Mass limit, the radius limit and index limit by a factor of 10
it will cause the program to take longer but will generate the star accurately
T_c range: 5e6 - 50e6
"""
T_0 = T_c = 5e6 # MS temperature paramater
X_0, Y_0, Z_0 = 0.7381, 0.2485, 0.0134   # metalicity  parameters change them here
mass_limit = 500 * 2e32 # the upper mass limit should be at least 10 solar masses and at most 1000
radius_limit = 20 * 6.957e8  # the upper radius limit should be at least 10 solar radii and at most 100
index_limit = 1e7 ## upper limit of index should be at least 1 billion

# --- Total Pressure, P (Equation 5) ---
"""
This function returns the sum of the contributions from non-relativistic
degenerate, ideal, and photon gases. As mentioned, this is a simple sum.

Required inputs:
    a = 
    rho = Density
    T = Temperature
    X, Y, Z = From metallicity
"""


def Total_P(rho, T, X, Y, Z):
    First_term = ((((3 * (pi ** 2)) ** (2 / 3.0)) * hbar ** 2) / (5 * m_e)) * ((rho / m_p) ** (5 / 3.0))
    Second_term = rho * ((k_b * T) / (mu(X, Y, Z) * m_p))
    Third_term = (1 / 3.0) * a * (T ** 4)

    return First_term + Second_term + Third_term


# --- Mean particle weight, mu (Equation 6) ---
"""
This function returns mu, given X, Y, Z
"""


def mu(X, Y, Z):
    return (2 * X + 0.75 * Y + 0.5 * Z) ** (-1)


# --- dP/d(rho) (Equation 7a) ---
"""
This function returns dP/d(rho).

Required inputs:
    rho = Density
    T = Temperature
    X, Y, Z = From metallicity
"""


def dP_drho(rho, T, X, Y, Z):
    First_term = ((((3 * (pi ** 2)) ** (2 / 3.0) * hbar ** 2)) / (3 * m_e * m_p)) * ((rho / m_p) ** (2 / 3.0))
    Second_term = (k_b * T) / (mu(X, Y, Z) * m_p)

    return First_term + Second_term


# --- dP/dT (Equation 7b)---
"""
This function returns dP/dT.

Required inputs:
    rho = Density
    T = Temperature
    X, Y, Z = From metallicity
"""


def dP_dT(rho, T, X, Y, Z):
    First_term = rho * (k_b) / (mu(X, Y, Z) * m_p)
    Second_term = (4 / 3.0) * a * (T ** 3)

    return First_term + Second_term

#--- Specific Energy Generation Rates (Equations 8a, 8b, 9) ---
"""
The functions below will output the specific energy generation rates in W/kg.
"""

#- Proton-Proton Chain (Equation 8a)-
def epsilon_pp(rho, T, X):
    return (1.07e-7) * (rho/1.0e5) * (X**2) * ((T/1.0e6)**4)

#- CNO (Equation 8b) -
def epsilon_CNO(rho, T, X):
    X_CNO = 0.03 * X
    return (8.24e-26) * (rho/1e5) * X * X_CNO * (T/1e6)**(19.9)

#- Total Specific Energy Generation Rate (Equation 9) -
def epsilon(rho, T, X):
    return (epsilon_pp(rho, T, X) + epsilon_CNO(rho, T, X))


# --- Rosseland Mean Opacities (Equations 10, 11, 12) ---
"""
This function returns the Rosseland mean opacities for different processes.
"""


def Kappa_es(X):
    return (0.02 * (1 + X))  # m**2/kg


def Kappa_ff(Z, rho, T):
    rho_3 = rho / (1.0e3)
    return ((1.0e24) * (Z + 0.0001) * (rho_3) ** (0.7) * T ** (-3.5))  # m**2/kg


def Kappa_H_minus(Z, rho, T):
    rho_3 = rho / (1.0e3)
    return ((2.5e-32) * (Z / 0.02) * (rho_3) ** (0.5) * T ** 9)  # m**2/kg


# --- Radiative Opacity (Equation 14) ---
"""
This function returns the radiative opacity
"""


def Kappa_Rad(rho, T, X, Z):
    # rho_3 = rho/(1.0e3)

    Kes = abs(Kappa_es(X))
    Kff = abs(Kappa_ff(Z, rho, T))

    if Kes > Kff:
        Kmax = Kes
    else:
        Kmax = Kff

    if Z == 0:
        return Kmax

    else:
        return 1 / ((1 / Kappa_H_minus(Z, rho, T)) + (1 / Kmax))\



# ---------- Differential Equations of Solve (Equation 2) ----------

"""
This section of the code is dedicated to Equation(s) 2, the differential
equations that we'll need to solve.
"""


# ----- Density Equations (2a) -----

def drho_dr(M, r, rho, dPdT, dTdr, dPdrho):
    return -1 * ((G * M * rho / (r ** 2)) + (dPdT * dTdr)) / (dPdrho)


# ----- Temperature Equations (2b) -----

# --- Radiative transfer dT/dr ---
"""
This function returns the radiative dT/dr, given r and T as inputs.
"""


def T_RT(L, r, T, kappa, rho):
    return (3 * kappa * rho * L) / (16 * pi * a * c * (T ** 3) * (r ** 2))


# --- Convective transfer dT/dr ---
"""
This function returns the convective dT/dr, given r and T as inputs.
"""


def T_CT(M, P, r, T, rho):
    return (1 - (1 / gamma)) * ((T * G * M * rho) / (P * (r ** 2)))


# --- dT/dr equation ---
"""
This function returns the "flatter" of T_RT and T_CT functions.
"""


def dT_dr(L, r, T, kappa, rho, M, P):
    TRT = abs(T_RT(L, r, T, kappa, rho))
    TCT = abs(T_CT(M, P, r, T, rho))

    if TRT < TCT:
        dT = TRT * (-1)
    else:
        dT = TCT * (-1)
    return dT


# --- Mass Equations (Equation 2c) ---

def dM_dr(r, rho):
    return 4 * pi * (r ** 2) * rho


# --- Luminosity equation (Equation 2d) ---
def dL_dr(r, rho, epsilon):
    return 4 * pi * (r ** 2) * rho * epsilon


# --- Optical depth equation (Equation 2e) ---
def dtau_dr(kappa, rho):
    return kappa * rho


# ----- Find Nearest -----
"""
Given an array, this function will find the closest value (specificed as "value"), and return its
position in the array. To find the value in the array, you MUST call array[returned value].
"""


def find_nearest(array, value):
    array = np.array(array)
    idx = (abs(array - value)).argmin()
    return array[idx]


# ----- Find Radius of Star -----

def find_radius_idx(tau, tau_inf):
    array = np.array(tau)
    idx = (abs(array - tau_inf)).argmin()
    return idx


# ----- Initial Conditions for Stellar Centre (Equation 15) -----
"""
This function will compute the initial conditions for the centre of a star from Equation 15.
NOTE 1: r_0 must be much, much smaller than any characteristic length (i.e. the assumed star
radius).
NOTE 2: We will need the hydrogen fraction for the specific energy production rate. Input the 
metallicity array, and this function will siphon off the X we need.
"""


def Initial_Conditions(rho_c, T_c, r_0, metallicity_array):
    X = metallicity_array[0]
    rho = rho_c
    T = T_c
    M = (4 * pi / 3) * ((r_0) ** 3) * rho_c
    L = M * epsilon(rho_c, T_c, X)

    return rho, T, M, L


# ----- Opacity Proxy (Equation 16) -----
"""
We know that the surface of the star is defined via opacity, specifically when 
tau(inf) - tau)R_star) =2/3. We will use Equation 16 as the opacity proxy: we want to use
tau(inf) - tau = kappa*rho**2/|dp/dr| = 2/3.

We wish to integrate until the opacity proxy is much less than 1.

tau(inf) is the tau at the largest radii considered.
tau(R_star) is then given by tau(inf)-(2/3)
"""


def Opacity_Proxy(drho_dr, kappa, rho):
    return kappa * (rho ** 2) / (abs(drho_dr))


# ----- Trial Solution Indicator (Equation 17) -----
"""
Refer to notes for the description of this function.
"""

# def f_error(L_star, R_star, T_star):

#    const =  (4 * pi * sigma * (R_star**2) * (T_star**4)
#    numerator = L_star - const
#    denominator = sqrt(const * L_star)

#    return numerator/denominator

def find_k1(rho, T, X, Y, Z, x, M, L, h):
    # k1

    ### evaluating the DE's to find the coefficients
    dPdT = dP_dT(rho, T, X, Y, Z)
    dPdrho = dP_drho(rho, T, X, Y, Z)
    kappa = Kappa_Rad(rho, T, X, Z)
    P = Total_P(rho, T, X, Y, Z)
    dTdr = dT_dr(L, x, T, kappa, rho, M, P)
    drhodr = drho_dr(M, x, rho, dPdT, dTdr, dPdrho)
    dMdr = dM_dr(x, rho)
    eps = epsilon(rho, T, X)
    dLdr = dL_dr(x, rho, eps)
    dtaudr = dtau_dr(kappa, rho)

    k_rho = h * drhodr
    k_T = h * dTdr
    k_M = h * dMdr
    k_L = h * dLdr
    k_tau = h * dtaudr

    return [k_rho, k_T, k_M, k_L, k_tau]


def find_k2(rho, T, X, Y, Z, x, M, L, h, k1):
    # k2
    k_rho1 = k1[0]
    k_T1 = k1[1]
    k_M1 = k1[2]
    k_L1 = k1[3]

    dx = h / 4.0
    drho = k_rho1 / 4.0
    dT = k_T1 / 4.0
    dM = k_M1 / 4.0
    dL = k_L1 / 4.0

    dPdT = dP_dT(rho + drho, T + dT, X, Y, Z)
    dPdrho = dP_drho(rho + drho, T + dT, X, Y, Z)
    kappa = Kappa_Rad(rho + drho, T + dT, X, Z)
    P = Total_P(rho + drho, T + dT, X, Y, Z)
    dTdr = dT_dr(L + dL, x + dx, T + dT, kappa, rho + drho, M + dM, P)
    drhodr = drho_dr(M + dM, x + dx, rho + drho, dPdT, dTdr, dPdrho)
    dMdr = dM_dr(x + dx, rho + drho)
    eps = epsilon(rho + drho, T + dT, X)
    dLdr = dL_dr(x + dx, rho + drho, eps)
    dtaudr = dtau_dr(kappa, rho + drho)

    k_rho = h * drhodr
    k_T = h * dTdr
    k_M = h * dMdr
    k_L = h * dLdr
    k_tau = h * dtaudr

    return [k_rho, k_T, k_M, k_L, k_tau]


def find_k3(rho, T, X, Y, Z, x, M, L, h, k1, k2):
    # k3

    k_rho1 = k1[0]
    k_T1 = k1[1]
    k_M1 = k1[2]
    k_L1 = k1[3]

    k_rho2 = k2[0]
    k_T2 = k2[1]
    k_M2 = k2[2]
    k_L2 = k2[3]

    dx = 3 * h / 8.0
    drho = (3 * (k_rho1) / 32.0) + (9 * k_rho2 / 32.0)
    dT = (3 * (k_T1) / 32.0) + (9 * k_T2 / 32.0)
    dM = (3 * (k_M1) / 32.0) + (9 * k_M2 / 32.0)
    dL = (3 * (k_L1) / 32.0) + (9 * k_L2 / 32.0)

    dPdT = dP_dT(rho + drho, T + dT, X, Y, Z)
    dPdrho = dP_drho(rho + drho, T + dT, X, Y, Z)
    kappa = Kappa_Rad(rho + drho, T + dT, X, Z)
    P = Total_P(rho + drho, T + dT, X, Y, Z)
    dTdr = dT_dr(L + dL, x + dx, T + dT, kappa, rho + drho, M + dM, P)
    drhodr = drho_dr(M + dM, x + dx, rho + drho, dPdT, dTdr, dPdrho)
    dMdr = dM_dr(x + dx, rho + drho)
    eps = epsilon(rho + drho, T + dT, X)
    dLdr = dL_dr(x + dx, rho + drho, eps)
    dtaudr = dtau_dr(kappa, rho + drho)

    k_rho = h * drhodr
    k_T = h * dTdr
    k_M = h * dMdr
    k_L = h * dLdr
    k_tau = h * dtaudr

    return [k_rho, k_T, k_M, k_L, k_tau]


def find_k4(rho, T, X, Y, Z, x, M, L, h, k1, k2, k3):
    # k4

    k_rho1 = k1[0]
    k_T1 = k1[1]
    k_M1 = k1[2]
    k_L1 = k1[3]

    k_rho2 = k2[0]
    k_T2 = k2[1]
    k_M2 = k2[2]
    k_L2 = k2[3]

    k_rho3 = k3[0]
    k_T3 = k3[1]
    k_M3 = k3[2]
    k_L3 = k3[3]

    dx = 12 * h / 13.0
    drho = (1932 * (k_rho1) / 2197.0) - (7200 * k_rho2 / 2197.0) + (7296 * k_rho3 / 2197.0)
    dT = (1932 * (k_T1) / 2197.0) - (7200 * k_T2 / 2197.0) + (7296 * k_T3 / 2197.0)
    dM = (1932 * (k_M1) / 2197.0) - (7200 * k_M2 / 2197.0) + (7296 * k_M3 / 2197.0)
    dL = (1932 * (k_L1) / 2197.0) - (7200 * k_L2 / 2197.0) + (7296 * k_L3 / 2197.0)

    dPdT = dP_dT(rho + drho, T + dT, X, Y, Z)
    dPdrho = dP_drho(rho + drho, T + dT, X, Y, Z)
    kappa = Kappa_Rad(rho + drho, T + dT, X, Z)
    P = Total_P(rho + drho, T + dT, X, Y, Z)
    dTdr = dT_dr(L + dL, x + dx, T + dT, kappa, rho + drho, M + dM, P)
    drhodr = drho_dr(M + dM, x + dx, rho + drho, dPdT, dTdr, dPdrho)
    dMdr = dM_dr(x + dx, rho + drho)
    eps = epsilon(rho + drho, T + dT, X)
    dLdr = dL_dr(x + dx, rho + drho, eps)
    dtaudr = dtau_dr(kappa, rho + drho)

    k_rho = h * drhodr
    k_T = h * dTdr
    k_M = h * dMdr
    k_L = h * dLdr
    k_tau = h * dtaudr

    return [k_rho, k_T, k_M, k_L, k_tau]


def find_k5(rho, T, X, Y, Z, x, M, L, h, k1, k2, k3, k4):
    # k5
    k_rho1 = k1[0]
    k_T1 = k1[1]
    k_M1 = k1[2]
    k_L1 = k1[3]

    k_rho2 = k2[0]
    k_T2 = k2[1]
    k_M2 = k2[2]
    k_L2 = k2[3]

    k_rho3 = k3[0]
    k_T3 = k3[1]
    k_M3 = k3[2]
    k_L3 = k3[3]

    k_rho4 = k4[0]
    k_T4 = k4[1]
    k_M4 = k4[2]
    k_L4 = k4[3]

    dx = h
    drho = (439 * (k_rho1) / 216.0) - (8 * k_rho2) + (3680 * k_rho3 / 513.0) - (845 * k_rho4 / 4104.0)
    dT = (439 * (k_T1) / 216.0) - (8 * k_T2) + (3680 * k_T3 / 513.0) - (845 * k_T4 / 4104.0)
    dM = (439 * (k_M1) / 216.0) - (8 * k_M2) + (3680 * k_M3 / 513.0) - (845 * k_M4 / 4104.0)
    dL = (439 * (k_L1) / 216.0) - (8 * k_L2) + (3680 * k_L3 / 513.0) - (845 * k_L4 / 4104.0)

    dPdT = dP_dT(rho + drho, T + dT, X, Y, Z)
    dPdrho = dP_drho(rho + drho, T + dT, X, Y, Z)
    kappa = Kappa_Rad(rho + drho, T + dT, X, Z)
    P = Total_P(rho + drho, T + dT, X, Y, Z)
    dTdr = dT_dr(L + dL, x + dx, T + dT, kappa, rho + drho, M + dM, P)
    drhodr = drho_dr(M + dM, x + dx, rho + drho, dPdT, dTdr, dPdrho)
    dMdr = dM_dr(x + dx, rho + drho)
    eps = epsilon(rho + drho, T + dT, X)
    dLdr = dL_dr(x + dx, rho + drho, eps)
    dtaudr = dtau_dr(kappa, rho + drho)

    k_rho = h * drhodr
    k_T = h * dTdr
    k_M = h * dMdr
    k_L = h * dLdr
    k_tau = h * dtaudr

    return [k_rho, k_T, k_M, k_L, k_tau]


def find_k6(rho, T, X, Y, Z, x, M, L, h, k1, k2, k3, k4, k5):
    # k6
    k_rho1 = k1[0]
    k_T1 = k1[1]
    k_M1 = k1[2]
    k_L1 = k1[3]

    k_rho2 = k2[0]
    k_T2 = k2[1]
    k_M2 = k2[2]
    k_L2 = k2[3]

    k_rho3 = k3[0]
    k_T3 = k3[1]
    k_M3 = k3[2]
    k_L3 = k3[3]

    k_rho4 = k4[0]
    k_T4 = k4[1]
    k_M4 = k4[2]
    k_L4 = k4[3]

    k_rho5 = k5[0]
    k_T5 = k5[1]
    k_M5 = k5[2]
    k_L5 = k5[3]

    dx = h / 2.0
    drho = (-8 * (k_rho1) / 27.0) + (2 * k_rho2) - (3544 * k_rho3 / 2565.0) + (1859 * k_rho4 / 4104.0) - (
                11 * k_rho5 / 40.0)
    dT = (-8 * (k_T1) / 27.0) + (2 * k_T2) - (3544 * k_T3 / 2565.0) + (1859 * k_T4 / 4104.0) - (11 * k_T5 / 40.0)
    dM = (-8 * (k_M1) / 27.0) + (2 * k_M2) - (3544 * k_M3 / 2565.0) + (1859 * k_M4 / 4104.0) - (11 * k_M5 / 40.0)
    dL = (-8 * (k_L1) / 27.0) + (2 * k_L2) - (3544 * k_L3 / 2565.0) + (1859 * k_L4 / 4104.0) - (11 * k_L5 / 40.0)

    dPdT = dP_dT(rho + drho, T + dT, X, Y, Z)
    dPdrho = dP_drho(rho + drho, T + dT, X, Y, Z)
    kappa = Kappa_Rad(rho + drho, T + dT, X, Z)
    P = Total_P(rho + drho, T + dT, X, Y, Z)
    dTdr = dT_dr(L + dL, x + dx, T + dT, kappa, rho + drho, M + dM, P)
    drhodr = drho_dr(M + dM, x + dx, rho + drho, dPdT, dTdr, dPdrho)
    dMdr = dM_dr(x + dx, rho + drho)
    eps = epsilon(rho + drho, T + dT, X)
    dLdr = dL_dr(x + dx, rho + drho, eps)
    dtaudr = dtau_dr(kappa, rho + drho)

    k_rho = h * drhodr
    k_T = h * dTdr
    k_M = h * dMdr
    k_L = h * dLdr
    k_tau = h * dtaudr

    return [k_rho, k_T, k_M, k_L, k_tau]


def RK45(i, rho, T, X, Y, Z, x, M, L, h, r, tau):
    hmin = 1e4
    hmax = 1e6
    htol = 1e-8
    step = h
    j = 0

    while j <= 1:
        k1 = find_k1(rho[i], T[i], X, Y, Z, x, M[i], L[i], step)
        k2 = find_k2(rho[i], T[i], X, Y, Z, x, M[i], L[i], step, k1)
        k3 = find_k3(rho[i], T[i], X, Y, Z, x, M[i], L[i], step, k1, k2)
        k4 = find_k4(rho[i], T[i], X, Y, Z, x, M[i], L[i], step, k1, k2, k3)
        k5 = find_k5(rho[i], T[i], X, Y, Z, x, M[i], L[i], step, k1, k2, k3, k4)
        k6 = find_k6(rho[i], T[i], X, Y, Z, x, M[i], L[i], step, k1, k2, k3, k4, k5)

        ### calculate 4th order runge kutta
        #### change this!!!

        rho4 = rho[i] + (25 * k1[0] / 216.0 + 1408 * k3[0] / 2565.0 + 2197 * k4[0] / 4104.0 - 1 * k5[0] / 5.0)
        T4 = T[i] + (25 * k1[1] / 216.0 + 1408 * k3[1] / 2565.0 + 2197 * k4[1] / 4104.0 - 1 * k5[1] / 5.0)
        M4 = M[i] + (25 * k1[2] / 216.0 + 1408 * k3[2] / 2565.0 + 2197 * k4[2] / 4104.0 - 1 * k5[2] / 5.0)
        L4 = L[i] + (25 * k1[3] / 216.0 + 1408 * k3[3] / 2565.0 + 2197 * k4[3] / 4104.0 - 1 * k5[3] / 5.0)
        tau4 = tau[i] + (25 * k1[4] / 216.0 + 1408 * k3[4] / 2565.0 + 2197 * k4[4] / 4104.0 - 1 * k5[4] / 5.0)

        ### calculate 5th order runge kutta
        rho5 = rho[i] + (
                    16 * k1[0] / 135.0 + 6656 * k3[0] / 12825.0 + 28561 * k4[0] / 56430.0 - 9 * k5[0] / 50.0 + 2 * k6[
                0] / 55.0)
        T5 = T[i] + (16 * k1[1] / 135.0 + 6656 * k3[1] / 12825.0 + 28561 * k4[1] / 56430.0 - 9 * k5[1] / 50.0 + 2 * k6[
            1] / 55.0)
        M5 = M[i] + (16 * k1[2] / 135.0 + 6656 * k3[2] / 12825.0 + 28561 * k4[2] / 56430.0 - 9 * k5[2] / 50.0 + 2 * k6[
            2] / 55.0)
        L5 = L[i] + (16 * k1[3] / 135.0 + 6656 * k3[3] / 12825.0 + 28561 * k4[3] / 56430.0 - 9 * k5[3] / 50.0 + 2 * k6[
            3] / 55.0)
        tau5 = tau[i] + (
                    16 * k1[4] / 135.0 + 6656 * k3[4] / 12825.0 + 28561 * k4[4] / 56430.0 - 9 * k5[4] / 50.0 + 2 * k6[
                4] / 55.0)

        diff = T5 - T4

        if j == 0:

            if diff == 0:
                s = 1
                j += 1
            elif isnan(htol / (2 * abs(diff)) ** (1 / 4.)):
                s = 0.5
            else:
                s = (htol / (2 * abs(diff)) ** (1 / 4.))

            step *= s

            if step > hmax:
                step = hmax
            elif step < hmin:
                step = hmin
        j += 1

    rho.append(rho5)
    T.append(T5)
    M.append(M5)
    L.append(L5)
    tau.append(tau5)
    r.append(x)

def RK5_fixed(i,rho,T,X,Y,Z,x,M,L,h,r,tau):
  k1 = find_k1(rho[i],T[i],X,Y,Z,x,M[i],L[i],h)
  k2 = find_k2(rho[i],T[i],X,Y,Z,x,M[i],L[i],h,k1)
  k3 = find_k3(rho[i],T[i],X,Y,Z,x,M[i],L[i],h,k1,k2)
  k4 = find_k4(rho[i],T[i],X,Y,Z,x,M[i],L[i],h,k1,k2,k3)
  k5 = find_k5(rho[i],T[i],X,Y,Z,x,M[i],L[i],h,k1,k2,k3,k4)
  k6 = find_k6(rho[i],T[i],X,Y,Z,x,M[i],L[i],h,k1,k2,k3,k4,k5)

  rho5 = rho[i] + (16*k1[0]/135.0 + 6656*k3[0]/12825.0 + 28561*k4[0]/56430.0 - 9*k5[0]/50.0 + 2*k6[0]/55.0)
  T5 = T[i] + (16*k1[1]/135.0 + 6656*k3[1]/12825.0 + 28561*k4[1]/56430.0 - 9*k5[1]/50.0 + 2*k6[1]/55.0)
  M5 = M[i] + (16*k1[2]/135.0 + 6656*k3[2]/12825.0 + 28561*k4[2]/56430.0 - 9*k5[2]/50.0 + 2*k6[2]/55.0)
  L5 = L[i] + (16*k1[3]/135.0 + 6656*k3[3]/12825.0 + 28561*k4[3]/56430.0 - 9*k5[3]/50.0 + 2*k6[3]/55.0)
  tau5 = tau[i] + (16*k1[4]/135.0 + 6656*k3[4]/12825.0 + 28561*k4[4]/56430.0 - 9*k5[4]/50.0 + 2*k6[4]/55.0)

  rho.append(rho5)
  T.append(T5)
  M.append(M5)
  L.append(L5)
  tau.append(tau5)
  r.append(x)


def f(rho_c):
    i = 0

    h = 1e4

    rho = [rho_c]
    T = [T_c]
    M = [M_0]
    L = [L_0]
    tau = [0]
    r = [0]
    x = 1e-15

    X, Y, Z = X_0, Y_0, Z_0

    tau_proxy = 1

    while (tau_proxy > 0.01) and (M[i - 1] < mass_limit) and (r[i - 1] < radius_limit):
        if i == 0:
            rho0, T0, Mi, Li = Initial_Conditions(rho_c, T_c, x, [X, Y, Z])
            print("rho_c =", rho_c)
            rho.append(rho_c)
            T.append(T_c)
            M.append(Mi)
            L.append(Li)
            kappa = Kappa_Rad(rho_c, T_c, X, Z)
            tau.append(0)
            r.append(x)
            i += 1
            x += h
        else:
            RK45(i, rho, T, X, Y, Z, x, M, L, h, r, tau)
            dPdT = dP_dT(rho[i], T[i], X, Y, Z)
            dPdrho = dP_drho(rho[i], T[i], X, Y, Z)
            kappa = Kappa_Rad(rho[i], T[i], X, Z)
            P = Total_P(rho[i], T[i], X, Y, Z)
            dTdr = dT_dr(L[i], x, T[i], kappa, rho[i], M[i], P)
            drhodr = drho_dr(M[i], x, rho[i], dPdT, dTdr, dPdrho)
            tau_proxy = Opacity_Proxy(drhodr, kappa, rho[i])
            if M[i] > 2e28 and i % 100000 == 0:
                if DEBUG_FLAG == 1:
                    print("i = " + str(i) + ", M = " + str(M[i] / 2e30) + " solar masses" + ", r = " + str(
                        r[i] / 6.955e8) + " solar radii")
            if M[i] > 2e32:
                if DEBUG_FLAG == 1:
                    print("stopped at i = " + str(i) + ", M = " + str(M[i] / 2e30) + " solar masses" + ", r = " + str(
                        r[i] / 6.955e8) + " solar radii")
            if tau_proxy < 0.01:
                if DEBUG_FLAG == 1:
                    print("converged at i = " + str(i) + ", M = " + str(M[i] / 2e30) + " solar masses" + ", r = " + str(
                        r[i] / 6.955e8) + " solar radii")

            x += h
            i += 1

    tau_inf = tau[-1] - 2 / 3.0
    idx = find_radius_idx(tau, tau_inf)

    L_s = L[idx]
    R_s = r[idx]
    T_s = T[idx]

    term1 = L_s - 4 * pi * sigma * (R_s ** 2) * (T_s ** 4)
    term2 = (4 * pi * sigma * (R_s ** 2) * (T_s ** 4) * L_s) ** (0.5)

    f = term1 / term2
    # print(f)

    return f


def final_star(rho_c):
    X, Y, Z = X_0, Y_0, Z_0
    h = 10000
    x = 1e-15
    rho = [rho_c]
    T = [T_c]
    M = [M_0]
    L = [L_0]
    tau = [0]
    P = [Total_P(rho_c, T_c, X, Y, Z)]
    kappa_l = [Kappa_Rad(rho_c, T_c, X, Z)]
    eps = epsilon(rho_c, T_c, X)
    dLdr = [dL_dr(x, rho_c, eps)]
    i = 0
    r = [0]

    file_name = str("Final_Star_T_c_" + str(T_c) + "_Z_" + str(Z) + ".csv")
    #  file_loc = str("C:\\Final_Star\\") + file_name
    file_loc = file_name
    col_names = ["rho", "T", "L", "Mass", "Tau", "Kappa_l", "P", "dL_dr", "r", "idx"]

    tau_proxy = 1

    while (tau_proxy > 0.01) and (M[i - 1] < mass_limit) and (r[i - 1] < radius_limit) and (i < 1e8):
        if i == 0:
            rho0, T0, Mi, Li = Initial_Conditions(rho_c, T_c, x, [X, Y, Z])
            rho.append(rho_c)
            T.append(T_c)
            M.append(Mi)
            L.append(Li)
            kappa = Kappa_Rad(rho_c, T_c, X, Z)
            tau.append(0)
            Press = Total_P(rho_c, T_c, X, Y, Z)
            P.append(Press)
            kappa_l.append(kappa)
            eps = epsilon(rho_c, T_c, X)
            dLdr.append(dL_dr(x, rho_c, eps))
            r.append(x)
            i += 1
            x += h

        else:
            RK45(i, rho, T, X, Y, Z, x, M, L, h, r, tau)
            dPdT = dP_dT(rho[i], T[i], X, Y, Z)
            dPdrho = dP_drho(rho[i], T[i], X, Y, Z)
            kappa = Kappa_Rad(rho[i], T[i], X, Z)
            Press = Total_P(rho[i], T[i], X, Y, Z)
            dTdr = dT_dr(L[i], x, T[i], kappa, rho[i], M[i], Press)
            drhodr = drho_dr(M[i], x, rho[i], dPdT, dTdr, dPdrho)
            tau_proxy = Opacity_Proxy(drhodr, kappa, rho[i])
            P.append(Press)
            kappa_l.append(kappa)
            eps = epsilon(rho[i], T[i], X)
            dLdr.append(dL_dr(x, rho[i], eps))
            if M[i] > 2e28 and i % 100000 == 0:
                if DEBUG_FLAG == 1:
                    print("i = " + str(i) + ", M = " + str(M[i] / 2e30) + " solar masses" + ", r = " + str(
                        r[i] / 6.955e8) + " solar radii")

            if tau_proxy < 0.01:
                if DEBUG_FLAG == 1:
                    print("converged at i = " + str(i) + ", M = " + str(M[i] / 2e30) + " solar masses" + ", r = " + str(
                        r[i] / 6.955e8) + " solar radii")

            x += h
            i += 1

    tau_inf = tau[-1] - 2 / 3.0
    idx = find_radius_idx(tau, tau_inf)

    L_s = L[idx]
    R_s = r[idx]
    T_s = T[idx]

    print("Save all of these parameters in txt file:")
    print("Z: ", Z_0)
    print("Core Temperature: ", T_c)
    print("Surface Temperature: ", T_s)
    print("Surface Luminosity: ", L_s)

    ### returns rho list, T list, L list, M list, tau list, kappa list, P list, dL_dr list, R list, index of R_star

    Data = pd.DataFrame(
        {'rho': rho,
         'T': T,
         'L': L,
         'M': M,
         'tau': tau,
         'kappa_l': kappa_l,
         'P': P,
         'dL/dr': dLdr,
         'r': r,
         'idx': idx
         })

    Data.to_csv(file_loc)

    return rho, T, L, M, tau, kappa_l, P, dLdr, r, idx


print("Starting")
T_list = np.linspace(5e6, 50e6, 10)
Ls = []
Ts = []

for T_l in T_list:
    print("T_c = ", T_l)
    T_c = T_l
    start_time = time()
    rho_c_ideal = bisect(f, 0.3e3, 500e3, xtol=1e-3)
    rho, T, L, M, tau, kappa, P, dLdr, r, idx = final_star(rho_c_ideal)
    end_time = time()
    run_time = end_time - start_time
    print("Total Run Time: ", run_time / 60.0)

    Ls.append(L[idx])
    Ts.append(T[idx])

plt.loglog(np.array(Ts),np.array(Ls)/3.828e26,'o')
plt.show
print(Ts,Ls)