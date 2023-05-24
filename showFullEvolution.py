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
import include.simulations.useQuasi3D as sim

plt.rcParams.update({'font.size': 15 })
mpl.use('Agg')


# Definition of Constants
M_E = 9.109e-31                      # Electron rest mass in kg
EC = 1.60217662e-19                  # Electron charge in C
EP_0 = 8.854187817e-12               # Vacuum permittivity in C/(V m)
C = 299892458                        # Speed of light in vacuum in m/s

# Snapshot locations (12 total, in mm):
#x_s = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 30]
#x_s = [500, 750, 1000, 1250, 400, 1750, 2000, 3000, 4000, 5000, 7500, 10000]
x_s = [0, 5, 10, 25, 50, 75, 100, 40, 200, 300, 400, 500]

# Color Scheme
WB = False # Sequential
Viridis = True # Sequential + Perceptually Uniform
BuPu = False # Sequential
Jet = False

t0 = sim.getTime()

def returnXi(z):
    return z - C * t0 

def returnZ(xi):
    return xi + C * t0

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

def plot(x_f,y_f,xi_f,z_f,px_f,py_f,pz_f, w, sim_name,shape_name,noElec,iter):
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
    slices = len(x_s) # Number of Screens
    xs_norm = []
    for i in range(0,slices):
        xs_norm.append(x_s[i] * W_P * 10**(-3) / C)

# Generate arrays of coordinates at origin + each screen
    yslice = np.empty([slices, noElec])
    xislice = np.empty([slices, noElec])
    zslice = np.empty([slices, noElec])

# Project positions at distances in x_s
    for i in range(0,slices):
        # If x_s out of plasma, use ballistic trajectory
        if (abs(xs_norm[i]) > plasma_bnds[2]):
            for j in range(0,noElec):
                yslice[i, j], xislice[i, j], zslice[i, j] = getBallisticTraj(x_f[j], y_f[j], xi_f[j], z_f[j], px_f[j], py_f[j], pz_f[j], xs_norm[i])
        else:
            for j in range(0,noElec):
                yslice[i, j] = y_f[j]
                xislice[i, j] = xi_f[j]
                zslice[i, j] = z_f[j]

# Plot slices
# For bin size = 0.006 (lambda/10)
# Run 130 Limits: (27,52), (-6,6), Bins: (4167,2000)
#         (35,40), (-1,1), Bins: (833,333)
# For bin size = 0.03
# Run 130 Limits: (27,52), (-6,6), Bins: (833,400)
# Run 232 Limits: (435,475), (0,6), Bins: (1333,200)

    binsizez = 1250#833#2833#4167#1000#2666#1333
    binsizey = 400#400#2000#160#666#200

    xmin = 25 #35#27#400
    xmax = 52 #500
    
    vmax_ = 10

    if (WB):
        cmap = plt.cm.binary
    elif (Viridis):
        cmap = plt.cm.viridis
    elif (BuPu):
        cmap = plt.cm.BuPu
    elif (Jet):
        cmap = copy.copy(plt.get_cmap('jet'))
        cmap.set_under(color='white')
    else:
        cmap = plt.cm.gist_gray
    norm = mpl.colors.Normalize(vmin=1, vmax=400)

    # Creating first plot
    fig5, axs = plt.subplots(3, sharey=True, figsize=(8, 10), dpi=600)
    fig5.suptitle("Progression of " + shape_name + " EProbe")
    for i in range(0, 3):
        #axs[i].set_title("X = " + str(x_s[i]) + " mm")
        h = axs[i].hist2d(zslice[i,:], yslice[i,:], weights=w, bins=(binsizez,binsizey), cmap=cmap, vmin=1, vmax=vmax_)#, norm=norm)
        axs[i].set_ylim(-6,6)
        axs[i].set_xlim(xmin,xmax)
        if (WB):
            axs[i].set_facecolor('white')
        elif (Viridis):
            axs[i].set_facecolor('#30013b')
        #elif (Jet):
            #axs[i].set_facecolor('#000040')
        else:
            axs[i].set_facecolor('white')

    axs[2].set(xlabel = 'Z ($c/\omega_p$)', ylabel = 'Y ($c/\omega_p$)')
    cbar = plt.colorbar(h[3], ax=axs, orientation='horizontal')
    # Creating secondary xi axis
    secax = axs[0].secondary_xaxis('top', functions= (returnXi, returnZ))
    secax.set(xlabel= '$\\xi$ ($c/\omega_p$)')
    #cbar.set_label('Electron Density')

    # Creating second plot
    fig6, axs2 = plt.subplots(3, sharey=True, figsize=(8, 10), dpi=600)
    fig6.suptitle("Progression of " + shape_name + " EProbe")
    for i in range(0, 3):
        #axs2[i].set_title("X = " + str(x_s[i+3]) + " mm")
        h2 = axs2[i].hist2d(zslice[i+3,:], yslice[i+3,:], weights=w, bins=(binsizez,int(binsizey/4)), cmin=0, cmap=cmap, vmin=1, vmax=vmax_)# norm=norm)
        axs2[i].set_ylim(-1.5,1.5)
        axs2[i].set_xlim(xmin,xmax)
        if (WB):
            axs2[i].set_facecolor('white')
        elif (Viridis):
            axs2[i].set_facecolor('#30013b')
        #elif (Jet):
            #axs2[i].set_facecolor('#000040')
        else:
            axs2[i].set_facecolor('white')

    axs2[2].set(xlabel = 'Z ($c/\omega_p$)', ylabel = 'Y ($c/\omega_p$)')
    cbar2 = plt.colorbar(h2[3], ax=axs2, orientation='horizontal')
    #cbar2.set_label('Electron Density')
    secax2 = axs2[0].secondary_xaxis('top', functions= (returnXi, returnZ))
    secax2.set(xlabel= '$\\xi$ ($c/\omega_p$)')

    fig7, axs3 = plt.subplots(3, sharey=True, figsize=(8, 10), dpi=600)
    fig7.suptitle("Progression of " + shape_name + " EProbe")
    for i in range(0, 3):
        #axs3[i].set_title("X = " + str(x_s[i+6]) + " mm")
        h3 = axs3[i].hist2d(zslice[i+6,:], yslice[i+6,:], weights=w, bins=(binsizez,int(binsizey/2)), cmap=cmap, vmin=1, vmax=vmax_)#, norm=norm)
        axs3[i].set_ylim(-3,3)
        axs3[i].set_xlim(xmin,xmax)
        if (WB):
            axs3[i].set_facecolor('white')
        elif (Viridis):
            axs3[i].set_facecolor('#30013b')
        #elif (Jet):
            #axs3[i].set_facecolor('#000040')
        else:
            axs3[i].set_facecolor('white')

    axs3[2].set(xlabel = 'Z ($c/\omega_p$)', ylabel = 'Y ($c/\omega_p$)')
    cbar3 = plt.colorbar(h3[3], ax=axs3, orientation='horizontal')
    #cbar3.set_label('Electron Density')
    secax3 = axs3[0].secondary_xaxis('top', functions= (returnXi, returnZ))
    secax3.set(xlabel= '$\\xi$ ($c/\omega_p$)')

    fig8, axs4 = plt.subplots(3, sharey=True, figsize=(8, 10), dpi=600)
    fig8.suptitle("Progression of " + shape_name + " EProbe")
    for i in range(0, 3):
        #axs4[i].set_title("X = " + str(x_s[i+9]) + " mm")
        if (i < 2):
            h4 = axs4[i].hist2d(zslice[i+9,:], yslice[i+9,:], weights=w, bins=(binsizez,binsizey), cmap=cmap, vmin=1, vmax=vmax_)#, norm=norm)
            axs4[i].set_ylim(-6,6)
            axs4[i].set_xlim(xmin,xmax)
        elif (i == 2):
            h4 = axs4[i].hist2d(zslice[i+9,:], yslice[i+9,:], weights=w, bins=(binsizez,binsizey), cmap=cmap, vmin=1, vmax=vmax_)#, norm=norm)
            axs4[i].set_ylim(-6,6)
            axs4[i].set_xlim(xmin,xmax)
        if (WB):
            axs4[i].set_facecolor('white')
        elif (Viridis):
            axs4[i].set_facecolor('#30013b')
        #elif (Jet):
            #axs4[i].set_facecolor('#000080')
        else:
            axs4[i].set_facecolor('white')

    axs4[2].set(xlabel = 'Z ($c/\omega_p$)', ylabel = 'Y ($c/\omega_p$)')
    cbar4 = plt.colorbar(h4[3], ax=axs4, orientation='horizontal')
    cbar4.set_label('Electron Density')
    secax4 = axs4[0].secondary_xaxis('top', functions= (returnXi, returnZ))
    secax4.set(xlabel= '$\\xi$ ($c/\omega_p$)')
    
    ##############################
    # Creating fifth plot
    fig9, axs5 = plt.subplots(3, sharey=True, figsize=(8, 10), dpi=600)
    fig9.suptitle("Progression of EProbe")
    for i in range(0, 3):
        #axs2[i].set_title("X = " + str(x_s[i+3]) + " mm")
        h5 = axs5[i].hist2d(zslice[i+3,:], yslice[i+3,:], weights=w, bins=(binsizez,binsizey), cmap=cmap)#, norm=norm)
        axs5[i].set_ylim(-1.5,1.5)
        axs5[i].set_xlim(xmin,xmax)
        if (WB):
            axs5[i].set_facecolor('white')
        elif (Viridis):
            axs5[i].set_facecolor('#30013b')
        else:
            axs5[i].set_facecolor('white')

    axs5[2].set(xlabel = 'Z ($c/\omega_p$)', ylabel = 'Y ($c/\omega_p$)')
    cbar5 = plt.colorbar(h5[3], ax=axs5, orientation='horizontal')
    #cbar5.set_label('Electron Density')
    secax5 = axs5[0].secondary_xaxis('top', functions= (returnXi, returnZ))
    secax5.set(xlabel= '$\\xi$ ($c/\omega_p$)')
    ################################

    #fig5.savefig('prog1.png',dpi=600,transparent=False)
    #fig6.savefig('prog2.png',dpi=600,transparent=False)
    #fig7.savefig('prog3.png',dpi=600,transparent=False)
    #fig8.savefig('prog4.png',dpi=600,transparent=False)
    fig9.savefig('prog5.png',dpi=600,transparent=False)
