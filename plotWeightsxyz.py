# Script for showing full evolution of probe at hardcoded snapshot locations in and out of plasma
import os
import numpy as np
import matplotlib.colors as col
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pdb
import math
import copy
import time
import progressbar
import include.movieWriter as movieWriter
import multiprocessing as mp
import include.simulations.useQuasi3D as sim

plt.rcParams.update({'font.size': 12 })
#plt.rcParams['animation.ffmpeg_path'] = '/ffmpeg/bin'
mpl.use('Agg')

# Definition of Constants
M_E = 9.109e-31                      # Electron rest mass in kg
EC = 1.60217662e-19                  # Electron charge in C
EP_0 = 8.854187817e-12               # Vacuum permittivity in C/(V m)
C = 299892458                        # Speed of light in vacuum in m/s

WB = False # Sequential
Viridis = True # Sequential + Perceptually Uniform
BuPu = False # Sequential
Jet = False

t0 = sim.getTime()

propspeed = sim.getPropagationSpeed()

def returnXi(z):
    return z - t0*propspeed

def returnZ(xi):
    return xi + t0*propspeed



def plotweightsxiy(y_f,xi_f, w, rand):
    
    path = os.getcwd()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    new_path = os.path.join(path,f'animation-{timestr}-{rand}')
    os.mkdir(new_path)

    plasma_bnds = sim.getBoundCond()
    shape_name = shape_name.capitalize()

    ximin = 15 #36  #25#27#400
    ximax = 35 #50  #500
    
    ymin = -6
    ymax = 6

    bin_resolution = 0.1 #0.02 #c/w_p
    bin_edges_xi = np.arange(ximin, ximax, bin_resolution)
    bin_edges_y = np.arange(ymin, ymax, bin_resolution)
    
    cmin = 1       # Minimum density displayed
    vmin_ = cmin    # Minimum color value
    vmax_ = 100    # Maximum color value

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
    
    # Create figure
    fig, ax = plt.subplots(1, figsize=(8, 5), dpi=600)
    fig.suptitle("Xi-Y Weighting Map")
    plt.tight_layout(rect=[0, 0, 1, 0.9])
    
    h = ax.hist2d(xi_f[:], y_f[:], weights=w[:], bins=(bin_edges_xi,bin_edges_y), cmap=cmap, vmin=vmin_,vmax=vmax_,cmin=cmin)#, norm=norm)

    ax.set_ylim(ymin,ymax)
    ax.set_xlim(ximin,ximax)

    if (WB):
        ax.set_facecolor('white')
    #elif (Viridis):
    #    ax.set_facecolor('#30013b')
    else:
        ax.set_facecolor('white')

    ax.set(xlabel = '$\\xi$ ($c/\omega_p$)', ylabel = 'Y ($c/\omega_p$)')

    secax = ax.secondary_xaxis('top', functions= (returnZ, returnXi))
    secax.set(xlabel= 'Z ($c/\omega_p$)')
    
    cbar = plt.colorbar(h[3], ax=ax, orientation='horizontal')#, pad=0.3)
    #cbar.set_label('Electron Density')

    #Saving
    filename = str(os.path.join(new_path,f'weighting-xi-y.png'))
    fig.savefig(filename,dpi=600,transparent=False)
        
    ax.cla()
    fig.clf()
    plt.close(fig)

    return

    
