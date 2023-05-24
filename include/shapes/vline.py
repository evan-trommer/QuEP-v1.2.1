import sys
import math
import numpy as np
import pdb

# Initalizes probe as vertical line of electrons with length 2*s1

def initProbe(x_c,y_c,xi_c,t0,s1,s2,s3,density, *args):

    x_0, y_0, xi_0, z_0 = [],[],[],[]

    step = 2*s1 / density
    yin = y_c - s1

    for i in range(0,density):
        x_0.append(x_c)
        xi_0.append(xi_c)
        z_0.append(xi_c + t0)
        y_0.append(yin)

        yin += step

    return x_0, y_0, xi_0, z_0
