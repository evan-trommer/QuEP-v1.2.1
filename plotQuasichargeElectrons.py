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

def getFieldArrays():

    xiaxis_1, xi2, raxis_1, r2 = sim.axes()
    xiiter = len(xiaxis_1)
    riter = len(raxis_1)

    chargeElectrons = np.empty((riter,xiiter),dtype=float)

    for ir in progressbar.progressbar(range(riter), redirect_stout=True):
        for ixi in range(xiiter):
            #pdb.set_trace()
            chargeElectrons[ir, ixi] = sim.chargeElectrons(1, raxis_1[ir], 0, xiaxis_1[ixi], raxis_1[ir], mode=0)

    return xiaxis_1, raxis_1, chargeElectrons

def main():

    start_time = time.time()
    t0 = sim.getTime()

    xiaxis, raxis, chargeElectrons = getFieldArrays()

    # Save data for future plotting
    fname = "chargeElectrons-plot-data.npz"
    np.savez(fname,xiaxis, raxis, chargeElectrons)
    print(f"Data for plot saved to {fname}")

    zaxis = [xi + t0 for xi in xiaxis]

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.8, top=0.9)
    fig.suptitle("charge Electrons, M0 Only")

    ax.set(xlabel = 'Z ($c/\omega_p$)', ylabel = 'X ($c/\omega_p$)')

    chargeElectrons = ax.pcolormesh(zaxis, raxis, chargeElectrons, norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-0.1,vmax=0.1),cmap="RdBu_r")

    tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
    cbar_ax = fig.add_axes([0.83, 0.05, 0.03, 0.9])

    cbar = fig.colorbar(chargeElectrons, cax=cbar_ax, ticks=tick_locations, format=ticker.LogFormatterMathtext())

    cbar.set_label('charge Electrons Density (Units ?)')

    print((time.time() - start_time)/60, " min")

    plt.savefig("chargeElectrons-density.png",transparent=True)
    #plt.show()

main()
