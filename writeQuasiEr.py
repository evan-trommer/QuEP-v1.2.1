import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pdb
import time
import include.simulations.useQuasi3D as sim
import csv

def getFieldArrays():

    xiaxis_1, xi2, raxis_1, r2 = sim.axes()
    xiiter = len(xiaxis_1)
    riter = len(raxis_1)

    Er_full = np.empty((riter,xiiter),dtype=float)
    Er_m0 = np.empty((riter,xiiter),dtype=float)
    Er_m1 = np.empty((riter,xiiter),dtype=float)

    for ir in range(riter):
        print(ir)
        for ixi in range(xiiter):
            #pdb.set_trace()
            Er_full[ir, ixi] = sim.EField(2, raxis_1[ir], 0, xiaxis_1[ixi], raxis_1[ir], mode=-1)
            Er_m0[ir, ixi] = sim.EField(2, raxis_1[ir], 0, xiaxis_1[ixi], raxis_1[ir], mode=0)
            Er_m1[ir, ixi] = sim.EField(2, raxis_1[ir], 0, xiaxis_1[ixi], raxis_1[ir], mode=1)

    return xiaxis_1, raxis_1, Er_full, Er_m0, Er_m1

def main():

    start_time = time.time()
    t0 = sim.getTime()

    xiaxis, raxis, Er_full, Er_m0, Er_m1 = getFieldArrays()
    zaxis = [xi + t0 for xi in xiaxis]

    with open('Axes.csv', 'w', newline='') as csvfile:
        axeswriter = csv.writer(csvfile, dialect='excel')
        axeswriter.writerow(xiaxis)
        axeswriter.writerow(zaxis)
        axeswriter.writerow(raxis)

    with open('Er_full.csv', 'w', newline='') as csvfile2:
        erfullwriter = csv.writer(csvfile2, dialect='excel')
        erfullwriter.writerows(Er_full)

    with open('Er_m0.csv', 'w', newline='') as csvfile3:
        erm0writer = csv.writer(csvfile3, dialect='excel')
        erm0writer.writerows(Er_m0)

    with open('Er_m1.csv', 'w', newline='') as csvfile4:
        erm1writer = csv.writer(csvfile4, dialect='excel')
        erm1writer.writerows(Er_m1)

    input()

main()
