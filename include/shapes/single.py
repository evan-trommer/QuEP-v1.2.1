import sys
import math
import numpy as np
import pdb

# Initializes a single electron at center coordinates

def initProbe(x_c, y_c, xi_c, t0, *args):
    x_0, y_0, xi_0, z_0 = [],[],[],[]
    x_0.append(x_c)
    y_0.append(y_c)
    xi_0.append(xi_c)
    z_0.append(xi_c + t0)

    return x_0, y_0, xi_0, z_0
