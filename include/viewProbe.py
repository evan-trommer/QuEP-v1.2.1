import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
import pdb
import math

def plot(x_dat,y_dat,xi_dat,z_dat,sim_name,shape_name,s1,s2,noElec):

    fig5 = plt.figure()
    ax5 = plt.axes()
    ax5.set_xlabel("$\\xi$ ($c/\omega_p$)")
    ax5.set_ylabel("Y ($c/\omega_p$)")
    ax5.set_title("Initial Electron Probe Shape")

    for i in range(0,noElec):
        ax5.scatter(xi_dat[i,0], y_dat[i,0], c='C1')

    fig5.show()
    input()
