# Script for showing full evolution of probe at hardcoded snapshot locations in and out of plasma

import numpy as np
import matplotlib.colors as col
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pdb
import math
import copy
import csv

# Definition of Constants
M_E = 9.109e-31                      # Electron rest mass in kg
EC = 1.60217662e-19                  # Electron charge in C
EP_0 = 8.854187817e-12               # Vacuum permittivity in C/(V m)
C = 299892458                        # Speed of light in vacuum in m/s

# Snapshot locations (12 total, in mm):
#x_s = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 30]
#x_s = [500, 750, 1000, 1250, 1500, 1750, 2000, 3000, 4000, 5000, 7500, 10000]
x_s = [0, 5, 10, 25, 50, 75, 100, 150, 200, 300, 400, 500]

def Gamma(p):
    return math.sqrt(1.0 + p**2)

def Velocity(px,ptot):
# Returns relativistic velocity from momentum
    return px / Gamma(ptot)

def getBallisticTraj(x_0,y_0,xi_0,z_0,px,py,pz,x_s):
# Use ballistic matrix to find positions on screens
    dx = x_s - x_0
    y_f = y_0 + dx * (py/px)
    z_f = z_0 + dx * (pz/px)

# Find time traveled to get proper xi
    p = math.sqrt(px**2 + py**2 + pz**2)
    vx = Velocity(px, p)
    vy = Velocity(py, p)
    vz = Velocity(pz, p)
    vtot = math.sqrt(vx**2 + vy**2 + vz**2)
    dtot = math.sqrt((x_s - x_0)**2 + (y_f - y_0)**2 + (z_f - z_0)**2)
    t = dtot/vtot

    xi_f = xi_0 + dx * (pz/px) + t

    return y_f, xi_f, z_f

def plot(x_f,y_f,xi_f,z_f,px_f,py_f,pz_f,sim_name,shape_name,noElec,iter):
# Plot evolution of probe after leaving plasma
    if (sim_name.upper() == 'OSIRIS_CYLINSYMM'):
        import include.simulations.useOsiCylin as sim
    elif (sim_name.upper() == 'QUASI3D'):
        import include.simulations.useQuasi3D as sim
    else:
        print("Simulation name unrecognized. Quitting...")
        exit()

    W_P = sim.getPlasFreq()
    plasma_bnds = sim.getBoundCond()
    shape_name = shape_name.capitalize()

# Normalize screen distances
    slices = len(x_s)
    xs_norm = []
    for i in range(0,slices):
        xs_norm.append(x_s[i] * W_P * 10**(-3) / C)

# Generate arrays of coordinates at origin + each screen
    yslice = np.empty(noElec)
    xislice = np.empty(noElec)
    zslice = np.empty(noElec)

# Project positions at distances in x_s
    for j in range(0,noElec):
        yslice[j], xislice[j] , zslice[j] = getBallisticTraj(x_f[j], y_f[j], xi_f[j], z_f[j], px_f[j], py_f[j], pz_f[j], xs_norm[-1])

# Plot slices
# For bin size = 0.006 (lambda/10)
# Run 130 Limits: (27,52), (-6,6), Bins: (4167,2000)
#         (35,40), (-1,1), Bins: (833,333)
# For bin size = 0.03
# Run 130 Limits: (27,52), (-6,6), Bins: (833,400)
# Run 232 Limits: (435,475), (0,6), Bins: (1333,200)

    binsizez = 833#2833#4167#1000#2666#1333
    binsizey = 400#2000#160#666#200

    xmin = 27#35#27#400
    xmax = 52#500

    norm = mpl.colors.Normalize(vmin=1, vmax=1500)

    h2, xbin, ybin, himage  = plt.hist2d(zslice, yslice, bins=(binsizez,binsizey), vmin=1)#, norm=norm)

    with open('counts.csv', 'w', newline='') as csvfile:
        nwriter = csv.writer(csvfile, dialect='excel')
        nwriter.writerows(h2)

    with open('xbins.csv', 'w', newline='') as csvfile2:
        binwriter = csv.writer(csvfile2, dialect='excel')
        binwriter.writerow(xbin)

    with open('ybins.csv', 'w', newline='') as csvfile3:
        binwriter2 = csv.writer(csvfile3, dialect='excel')
        binwriter2.writerow(ybin)
