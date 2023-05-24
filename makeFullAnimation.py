# Script for showing full evolution of probe at hardcoded snapshot locations in and out of plasma
import os
import numpy as np
import matplotlib.colors as col
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from matplotlib import rc
import pdb
import math
import copy
import time
import progressbar
import multiprocessing as mp
import include.simulations.useQuasi3D as sim

plt.rcParams.update({'font.size': 12 })
#plt.rcParams['animation.ffmpeg_path'] = '/ffmpeg/bin'
plt.rcParams['figure.constrained_layout.use'] = True
mpl.use('Agg')

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

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

def Gamma(p):
    return np.sqrt(1.0 + p**2)

def Velocity(px,ptot):
# Returns relativistic velocity from momentum
    return px / Gamma(ptot)

def getBallisticTraj(x_0,y_0,z_0,px,py,pz,x_s):
# Use ballistic matrix to find positions on screens
    dx = x_s - x_0
    y_f = y_0 + dx * (py/px)
    z_f = z_0 + dx * (pz/px)

# Find time traveled to get proper xi
    p = np.sqrt(px**2 + py**2 + pz**2)
    vx = Velocity(px, p)
    vy = Velocity(py, p)
    vz = Velocity(pz, p)
    vtot = np.sqrt(vx**2 + vy**2 + vz**2)
    dtot = np.sqrt((x_s - x_0)**2 + (y_f - y_0)**2 + (z_f - z_0)**2)
    t = dtot/vtot

    #xi_f = xi_0 + dx * (pz/px) + t

    return y_f, z_f



def prepare(sim_name,shape_name,noObj,rand):
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

# Plot slices
# For bin size = 0.006 (lambda/10)
# Run 130 Limits: (27,52), (-6,6), Bins: (4167,2000)
#         (35,40), (-1,1), Bins: (833,333)
# For bin size = 0.03
# Run 130 Limits: (27,52), (-6,6), Bins: (833,400)
# Run 232 Limits: (435,475), (0,6), Bins: (1333,200)

    ######## PLOT PARAMETERS: ########
    
    # Choose boundaries of screens in mm
    xstart_mm = 0
    xend_mm = 31
    xstep_mm = 31

    #binsizez = 6500//4#6000#833#2833#4167#1000#2666#1333
    #binsizey = 1000//4#400#2000#160#666#200
    
    # For Quasi_ID = 000130, use (36,48) for front blowout regime, or (28,50) for whole wakefield structure
    # For Quasi_ID = 000067, use (24,37)
    ximin = -14.5
    ximax = -4.5
    
    zmin =  ximin+t0#28  #25#27#400
    zmax =  ximax+t0#50  #500
    
    ymin = -0.5
    ymax = 0.5

    bin_resolution = 0.02 #0.02 #c/w_p
    bin_edges_z = np.arange(zmin, zmax, bin_resolution)
    bin_edges_y = np.arange(ymin, ymax, bin_resolution)
    
    cmin = 1       # Minimum density displayed
    vmin_ = cmin    # Minimum color value
    vmax_ = 100    # Maximum color value

    fps = 4 # frames per second for movie

    ######## END PARAMETERS ########

    if (WB):
        cmap = plt.cm.binary
    elif (Viridis):
        cmap = plt.cm.plasma#viridis; 'YlOrRd'
    elif (BuPu):
        cmap = plt.cm.BuPu
    elif (Jet):
        cmap = copy.copy(plt.get_cmap('jet'))
        cmap.set_under(color='white')
    else:
        cmap = plt.cm.gist_gray
    norm = mpl.colors.Normalize(vmin=1, vmax=400)

    # Normalize screen distances
    screen_dists = list(range(xstart_mm,xend_mm+1,xstep_mm))#np.arange(xstart_mm,xend_mm+1,xstep_mm)#
    print(f"First screen location in mm is x = {screen_dists[0]}, and last screen is at x = {screen_dists[-1]}")
    slices = len(screen_dists) # Number of Screens
    xs_norm = []
    for i in range(0,slices):
        xs_norm.append(screen_dists[i] * W_P * 10**(-3) / C)
    print(f"First screen located at x = {xs_norm[0]} and last screen recorded at x = {xs_norm[-1]}")

    # Generate arrays of coordinates at origin + each screen
    yslice = np.empty([noObj])
    zslice = np.empty([noObj])

    # Get cwd and create path variable for frame output
    path = os.getcwd()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    new_path = os.path.join(path,f'probe-prop-{timestr}-{rand}')
    #os.mkdir(new_path)

    return plasma_bnds, slices, xs_norm, yslice, zslice, bin_edges_z, bin_edges_y, cmap, cmin, vmin_, vmax_, zmin, zmax, ymin, ymax, fps, new_path, screen_dists
    


def plotmp(i,x_f,y_f,z_f,px_f,py_f,pz_f, w, xden, plasma_bnds, xs_norm, yslice, zslice, bin_edges_z, bin_edges_y, cmap, cmin, vmin_, vmax_, zmin, zmax, ymin, ymax, new_path, screen_dists):
    # Create figure
    fig, ax = plt.subplots(1, figsize=(8, 5), dpi=600)
    #fig.suptitle("Progression of Electron Probe")
    #plt.tight_layout(rect=[0, 0, 1, 0.9])
        
    # Project positions at distances in x_s
    # If x_s out of plasma, use ballistic trajectory
    #print(xs_norm[0], xs_norm[-1])
    if (abs(xs_norm[i]) > plasma_bnds[2]):
        yslice, zslice = getBallisticTraj(x_f, y_f, z_f, px_f, py_f, pz_f, xs_norm[i])
    else:
        yslice = y_f
        zslice = z_f
    
    timestr = time.strftime("%Y%m%d-%H%M%S")
    
    h = ax.hist2d(zslice[:], yslice[:], weights=w[:], bins=(bin_edges_z,bin_edges_y), cmap=cmap, vmin=vmin_,vmax=vmax_,cmin=cmin)#, norm=norm)
    temptext = ax.text(zmin+2,ymax*0.75,f"x = {screen_dists[i]:03}mm", fontdict=None, horizontalalignment='left', fontsize=18, color="Black")
    
    ax.set_ylim(ymin,ymax)
    ax.set_xlim(zmin,zmax)
    if (WB):
        ax.set_facecolor('white')
    #elif (Viridis):
    #    ax.set_facecolor('#30013b')
    else:
        ax.set_facecolor('white')

    #ax.set(xlabel = 'Z ($c/\omega_p$)', ylabel = 'Y ($c/\omega_p$)')
    fig.suptitle(fr'Probe Propagation', size = 22)

    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    
    ax.set_xlabel(r'Z ($c/\omega_p$)', size = 20)
    ax.set_ylabel(r'Y ($c/\omega_p$)', size = 20)
    
    ax.set_xticks([38,46])
    ax.set_yticks([-0.5,-0.3,0.3,0.5])

    ax.xaxis.set_label_coords(0.5,-0.02)
    ax.yaxis.set_label_coords(-0.01,0.5)
    
    #secax = ax.secondary_xaxis('top', functions= (returnXi, returnZ))
    #secax.set(xlabel= '$\\xi$ ($c/\omega_p$)')
    
    cbar = plt.colorbar(h[3], ax=ax, orientation='horizontal')#, pad=0.2)
    #cbar.set_label('Electron Density')

    #Saving
    filenumber = "{:05.1f}".format(screen_dists[i]).replace(".","-")
    filename = str(os.path.join(new_path,f'progression-x-{filenumber}mm-5.png'))
    fig.savefig(f'progression-x-{filenumber}mm-{timestr}.png',dpi=600,transparent=False)
        
    ax.cla()
    fig.clf()
    plt.close(fig)
    
    #temptext.remove()
    #cbar.remove()
    ################################

    
