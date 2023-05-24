import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pdb
import statistics as stat
import math

def calculate(x_0,y_0,xi_0,z_0,x_dat,y_dat,z_dat,xi_dat,px_f,py_f,pz_f,sim_name,shape_name,x_s,s1,s2):

    y0 = np.array(y_0)
    pxf = np.array(px_f)
    pyf = np.array(py_f)
    xdat = np.array(x_dat)
    ydat  = np.array(y_dat)

    focal = abs(y0) * pxf / abs(pyf)

    avgFocal = np.mean(focal)
    stdFocal = np.std(focal)

    searchRange = np.where(xdat[0,:] > avgFocal - stdFocal or xdat[0,:] < avgFocal + stdFocal)
    print(searchRange)