# Script for generating 2D plots of electron trajectories with option for plotting force

import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import rc
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
import pdb

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

from numpy.core.fromnumeric import size
plt.rcParams.update({'font.size': 16})

plotYForce = True # Plot transverse force with trajectories, not useful for many trajectories
plotZForce = True # Plot force along WF propagation

#large_size = 12

#plt.rc('ytick', labelsize=large_size)
#plt.rc('axes', labelsize=large_size)

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def plot(x_dat,y_dat,z_dat,xi_dat,Fx_dat,Fy_dat,Fz_dat,px_dat,py_dat,pz_dat,sim_name,shape_name,s1,s2,noElec,fname):

    
# 2D: Z-X, constrained to blowout regime
    fig1 = plt.figure(1, figsize=(15,10),dpi=300)
    ax1 = plt.axes()
    ax1.set_xlabel("X ($c/\omega_p$)")
    ax1.set_ylabel("Z ($c/\omega_p$)")
    ax1.set_xlim(-1,2)
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.set_title("Electron Trajectories through Blowout Regime")

    for i in range(0, noElec):
        ax1.plot(x_dat[i,:], z_dat[i,:], 'k', label='Z-X Trajectory') # Want vertical axis as y

    if (plotZForce):
        ax1_f = ax1.twinx()
        ax1_f.set_ylabel("$F_z$ ($m_e c \omega_p$)")
        ax1_f.yaxis.label.set_color('C0')
        ax1_f.tick_params(axis='y', labelcolor='C0', colors='C0')

        for i in range(0, noElec):
            ax1_f.plot(x_dat[i,:], Fz_dat[i,:], 'C0', label='Z Force')

        fig1.legend(bbox_to_anchor=(0.88, 0.94), bbox_transform=plt.gcf().transFigure)

# 2D: Y-X
    
    #print(type(px_dat), type(px_dat[0]), len(x_dat[0]), len(px_dat[0]))
    
    def getBallisticTraj1(x_dat,y_dat,px_dat,py_dat,pz_dat, noElec):
    # Use ballistic matrix to find positions on screens
        screen_min = x_dat[0, -1]
        screen_max = 4000
        x_steps = 8000
        dx = (screen_max - screen_min)/x_steps
        xpoints = np.linspace(screen_min, screen_max, x_steps)
        
        x_f = np.empty((noElec, x_steps))
        y_f = np.empty((noElec, x_steps))
        
        for i in range(0, noElec): 
            x_f[i][0] = screen_min
            y_f[i][0] = y_dat[0, -1]
        
        for i in range(0, noElec):
            for k in range(1, len(xpoints)): 
                x_f[i][k] = xpoints[k]
                y_f[i][k] = y_f[i][k-1] + dx * (py_dat[i][-1]/px_dat[i][-1])
        
        return x_f, y_f

    fig2, ax2 = plt.subplots(1,figsize=(15,10),dpi=300)
    fig2.subplots_adjust(right=0.7)

    #print(f"X range: {x_dat[0,0]} to {x_dat[0,-1]}")
    #print(f"Y range: {y_dat[0,0]} to {y_dat[0,-1]}")
    
    x_coord, y_coord = getBallisticTraj1(x_dat,y_dat,px_dat,py_dat,pz_dat,noElec)    
    x_dat1 = np.empty(( noElec, (len(x_coord[0]) + len(x_dat[0])) ))
    y_dat1 = np.empty(( noElec, (len(x_coord[0]) + len(x_dat[0])) ))
    
    for i in range(0, noElec):
        x_dat1[i] = np.concatenate((x_dat[i],x_coord[i]))
        y_dat1[i] = np.concatenate((y_dat[i],y_coord[i]))
    
    for i in range(0, noElec):
        #y_dat[i,:] = [y/0.65 for y in y_dat[i,:]]
        #ax2.plot(x_dat[i,:], y_dat[i,:], 'k', label='Y-X Electron Trajectory')
        ax2.plot(x_dat1[i,:], y_dat1[i,:], 'k', label='Y-X Electron Trajectory')# Want vertical axis as y
        ax2.set_xlim(-1.5,1.5)
        ax2.set_ylim(-0.1, 0.65)
        #ax2.axhline(0, 0, 1, color = 'b', ls = '--', lw = 0.8)
    #ax2.axhline(y = 0, xmin=0, xmax=1, color = 'b')
    #ax2.axvline(x = 407, ymin = 0, ymax=1, color = 'b', ls = '--', lw = 0.8)
    ax2.axvline(x = -0.55, ymin=0, ymax=1, ls = '--', color = 'b')
    ax2.axvline(x = 0.6, ymin=0, ymax=1, ls = '--', color = 'b')
    ax2.set_xlabel("X ($c/\omega_p$)")
    ax2.set_ylabel("Y ($c/\omega_p$)")
    ax2.set_title("Electron Trajectory through Blowout Regime")

    if (plotYForce):
        Fy_ax = ax2.twinx()
        px_ax = ax2.twinx()
        py_ax = ax2.twinx()
        #pz_ax = ax2.twinx()

        px_ax.spines["right"].set_position(("axes",1.05))
        make_patch_spines_invisible(px_ax)
        px_ax.spines["right"].set_visible(True)
        py_ax.spines["right"].set_position(("axes",1.20))
        make_patch_spines_invisible(py_ax)
        py_ax.spines["right"].set_visible(True)
        #pz_ax.spines["right"].set_position(("axes",1.35))
        #make_patch_spines_invisible(pz_ax)
        #pz_ax.spines["right"].set_visible(True)
        Fy_ax.spines["right"].set_position(("axes",1.35))
        make_patch_spines_invisible(Fy_ax)
        Fy_ax.spines["right"].set_visible(True)

        for i in range(0, noElec):
            Fy_ax.plot(x_dat[i,:], Fy_dat[i,:], 'C0', label='Transverse Electric Force, $F_y$')
            px_ax.plot(x_dat[i,:], px_dat[i,:], 'C1', label='Momentum in X')
            py_ax.plot(x_dat[i,:], py_dat[i,:], 'C2', label='Momentum in Y')
            #pz_ax.plot(x_dat[i,:], pz_dat[i,:], 'C0', label='Momentum in Z')
        Fy_ax.set_ylabel("$F_y$ ($m_e c \omega_p$)")
        px_ax.set_ylabel("$p_x (m_e c)$")
        py_ax.set_ylabel("$p_y (m_e c)$")
        #pz_ax.set_ylabel("$p_z (m_e c)$")

        Fy_ax.yaxis.label.set_color('C0')
        px_ax.yaxis.label.set_color('C1')
        py_ax.yaxis.label.set_color('C2')
        #pz_ax.yaxis.label.set_color('C0')

        tkw = dict(size=4, width=1.5)
        #ax2.tick_params(axis='y', colors='k', **tkw)
        Fy_ax.tick_params(axis='y', colors='C0', **tkw)
        px_ax.tick_params(axis='y', colors='C1', **tkw)
        py_ax.tick_params(axis='y', colors='C2', **tkw)
        #pz_ax.tick_params(axis='y', colors='C0', **tkw)
        ax2.tick_params(axis='x', **tkw)

        ax2.grid()
        fig2.legend(bbox_to_anchor=(0.3, 0.8), bbox_transform=plt.gcf().transFigure)

    #fig1.tight_layout()
    #fig1.savefig(f"eProbe-xzTrajectories_{fname}.png",transparent=False)
    #fig2.tight_layout()
    fig2.savefig(f"eProbe-xyTrajectories_{fname}.png",transparent=False)
