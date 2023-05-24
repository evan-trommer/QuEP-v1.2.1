import sys
import math
import numpy as np
import pdb

# Initializes probe as a circle of electrons with radius s1

# Constants
PI = math.pi

def initProbe(x_c,y_c,xi_c,t0,s1,s2,s3,density):
    x_0, y_0, xi_0, z_0 = [],[],[],[]
    step = (2 * PI) / density
    theta = 0

    while (theta <= 2*PI):
        x_0.append(x_c)
        xi_0.append(xi_c + s1*math.cos(theta))
        y_0.append(y_c + s1*math.sin(theta))
        z_0.append(xi_c + s1*math.cos(theta) + t0)
        theta += step

    return x_0, y_0, xi_0, z_0
