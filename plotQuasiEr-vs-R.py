import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pdb
import time
import progressbar
import include.simulations.useQuasi3D as sim

xi_position = [-8.3] # Xi Coordinates for Er vs R slices
xilength = len(xi_position)

def getFieldArrays():

    xiaxis_1, xi2, raxis_1, r2 = sim.axes()
    
    riter = len(raxis_1)

    Er_m0 = np.empty((riter,xilength),dtype=float)

    for ir in progressbar.progressbar(range(riter), redirect_stout=True):
        for ixi in range(xilength):
            Er_m0[ir,ixi] = sim.EField(2, raxis_1[ir], 0, xi_position[ixi], raxis_1[ir], mode=0)

    return xiaxis_1, raxis_1, Er_m0

def main():

    start_time = time.time()
    t0 = sim.getTime()

    xiaxis, raxis, Er_m0 = getFieldArrays()

    # Save data for future plotting
    fname = "Er-vs-R-plot-data.npz"
    np.savez(fname,xiaxis, raxis, Er_m0)
    print(f"Data for plot saved to {fname}")

    zaxis = [xi + t0 for xi in xiaxis]

    fit_max = 100
    k, r0 = np.polyfit(raxis[:fit_max], Er_m0[:fit_max,0],1)


    fig1, ax1 = plt.subplots(figsize=(10,8))

    fig1.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9)

    Er_fit = np.empty((len(raxis)),dtype=float)
    for ir in range(len(raxis)):
        Er_fit[ir] = k*raxis[ir] + r0

    for ixi in range(xilength):
        ax1.plot(raxis[:], Er_m0[:,ixi], label=f"Xi = {xi_position[ixi]}")
        ax1.plot(raxis[:fit_max], Er_m0[:fit_max,ixi], label=f"Xi = {xi_position[ixi]}",color="orange")
        ax1.plot(raxis[:fit_max], Er_fit[:fit_max], "--", label=f"Er fit, Er={k:.3f}*r+{r0:.3f}",color="green")
        ax1.axvline(x = 0.65, ymin = 0, ymax = 1, color = 'k', ls = '--')
    
    ax1.set(xlabel = 'R ($c/\omega_p$)', ylabel = '$E_r$ ($c/\omega_p$)')

    ax1.set_title('Radial Electric Field as a Function of R, M0 Only')

    ax1.legend()

    print((time.time() - start_time)/60, " min")

    fig1.savefig("Er-vs-R-M0-fields.png",dpi=600,transparent=True)

main()
