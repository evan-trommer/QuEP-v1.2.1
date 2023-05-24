import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pdb
import time
import progressbar
import include.simulations.useQuasi3D as sim
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def getFieldArrays():

    xiaxis_1, xi2, raxis_1, r2 = sim.axes()
    xiiter = len(xiaxis_1)
    riter = len(raxis_1)

    Ez = np.empty((riter,xiiter),dtype=float)

    for ir in progressbar.progressbar(range(riter), redirect_stout=True):
        for ixi in range(xiiter):
            #pdb.set_trace()
            Ez[ir, ixi] = sim.EField(1, raxis_1[ir], 0, xiaxis_1[ixi], raxis_1[ir], mode=0)

    return xiaxis_1, raxis_1, Ez

def main():

    start_time = time.time()
    t0 = sim.getTime()

    xiaxis, raxis, Ez = getFieldArrays()

    # Save data for future plotting
    fname = "Ez-plot-data.npz"
    np.savez(fname,xiaxis, raxis, Ez)
    print(f"Data for plot saved to {fname}")

    def getXi(zaxis): 
        xiaxis = zaxis -t0 
        return xiaxis 
    
    def getZ(xiaxis): 
        zaxis = xiaxis + t0
        return zaxis
    
    #zaxis = [xi + t0 for xi in xiaxis]

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.8, top=0.85)
    fig.suptitle("Transverse ($\\phi = 0$) Electric Field in Z, M0 Only")

    ax.set(xlabel = r'$\xi$ ($c/\omega_p$)', ylabel = 'X ($c/\omega_p$)')
    secax = ax.secondary_xaxis('top', functions=(getZ, getXi))
    secax.set_xlabel('Z ($c/\omega_p$)')

    Ez = ax.pcolormesh(xiaxis, raxis, Ez, vmin = -0.1, vmax=0.1, cmap="RdBu_r")#, norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-0.1,vmax=0.1),cmap="RdBu_r")
    ax.set_ylim(0,2.0)

    #ax.axhline(y = 0.05, xmin=0, xmax=1, color="k", ls="--")
    #ax.axhline(y = 0.1, xmin=0, xmax=1, color="k", ls="--")
    #ax.axhline(y = 0.25, xmin=0, xmax=1, color="k", ls="--")
    #ax.axhline(y = 0.5, xmin=0, xmax=1, color="k", ls="--")
    #ax.axhline(y = 0.65, xmin=0, xmax=1, color="k", ls="--")
    #ax.axhline(y = 1, xmin=0, xmax=1, color="k", ls="--")
    
    ax.axvline(x = -5.5, ymin = 0, ymax = 1, color='k')
    ax.axvline(x = -15.5, ymin = 0, ymax = 1, color='k')
    
    #tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
    cbar_ax = fig.add_axes([0.83, 0.05, 0.03, 0.9])

    cbar = fig.colorbar(Ez, cax=cbar_ax, extend = "both")#, ticks=tick_locations, format=ticker.LogFormatterMathtext())

    cbar.set_label('Electric Field ($m_e c \omega_p / e$)')

    print((time.time() - start_time)/60, " min")

    plt.savefig("Ez-fields.png",dpi=600, transparent=True)
    #plt.show()
    input()

main()
