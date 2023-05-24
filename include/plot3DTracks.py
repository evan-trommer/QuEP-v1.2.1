# Script for generating 3D plots of electron trajectories with option for plotting force

import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
import pdb

def plot(x_dat,y_dat,xi_dat,z_dat,x_f,y_f,xi_f,z_f,px_f,py_f,pz_f,sim_name,shape_name,s1,s2,noElec):

# 3D: X, Xi, Y
    shape_name = shape_name.capitalize()
    fig3 = plt.figure(3)
    ax3 = plt.axes(projection='3d')
    ax3.set_xlabel("X ($c/\omega_p$)")
    ax3.set_ylabel("$\\xi$ ($c/\omega_p$)")
    ax3.set_zlabel("Y ($c/\omega_p$)")
    ax3.set_title(shape_name + " Electron Probe Trajectories in $\\xi$")
    for i in range(0, noElec):
        ax3.plot(x_dat[i,:], xi_dat[i,:], y_dat[i,:], 'k') # Want vertical ax3is as ya

    #ax3.legend(bbox_to_anchor=(0.97, 0.97), bbox_transform=plt.gcf().transFigure)

# 3D: X, Z, Y
    fig4 = plt.figure(4)
    ax4 = plt.axes(projection='3d')
    ax4.set_xlabel("X ($c/\omega_p$)")
    ax4.set_ylabel("Z ($c/\omega_p$)")
    ax4.set_zlabel("Y ($c/\omega_p$)")
    ax4.set_title(shape_name + " Electron Probe Trajectories in $\\xi$")
    for i in range(0, noElec):
        ax4.plot(x_dat[i,:], z_dat[i,:], y_dat[i,:], 'k') # Want vertical ax3is as y
    #ax4.legend(bbox_to_anchor=(0.97, 0.97), bbox_transform=plt.gcf().transFigure)

    fig3.show()
    fig4.show()
    input()
