import sys
import math
import numpy as np
import pdb

# Initalizes probe as horizontal line of electrons with length 2*s2

def initProbe(x_c,y_c,xi_c,t0,s1,s2,s3,density):

    x_0, y_0, xi_0, z_0 = [],[],[],[]

    step = 2*s2 / density
    xin = xi_c - s2

    for i in range(0,density):
        x_0.append(x_c)
        y_0.append(y_c)
        xi_0.append(xin)
        z_0.append(xin + t0)

        xin += step

    return x_0, y_0, xi_0, z_0
