import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pdb
import statistics as stat
import math

# Initialize a vertical line of electrons through the blowout regime to find the Y focal length
# Spherical aberrations will cause variance in focal length

rb = 0.65
k = 0.19#0.493#0.475
vx = 1 # =0.1 when pz = 1000, =1 when pz = 0

def calculate(x_0,y_0,xi_0,z_0,x_f,y_f,xi_f,z_f,px_f,py_f,pz_f,sim_name,shape_name,x_s,s1,s2):

    dy, fAngle, fXp = [],[],[]

# With thin lens approximation, dy = y0 - 0
    for i in range(0, len(x_0)):
        dy.append(y_0[i])

# Find focal length for each data point by multiplying dy by the angle of the electron leaving the regime
    for i in range(0, len(x_0)):
        fAngle.append( abs(dy[i]) * px_f[i] / abs(py_f[i]))
        fXp.append(px_f[i] * vx / (2 * math.sqrt(rb**2 - y_0[i]**2) * k ))
        print("px_f = ", px_f[i], ", py_f = ", py_f[i])
        print("f from p_y = ", fAngle[i])
        print("f from x_p = ", fXp[i])
        #pdb.set_trace()
# Find average focal length and variance
    if (len(dy) > 1):
        focal_y = stat.mean(fAngle)
        std_y = math.sqrt(stat.variance(fAngle))
        print("Focal Y from p_y = " + str(focal_y) + " " + u"\u00B1 " + str(std_y))
    else:
        print("Focal Y from p_y= " + str(fAngle))
    input()