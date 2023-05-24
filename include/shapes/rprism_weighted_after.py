# IMPORTANT Note:
# This file is depricated and is only still here for future reference
# 3D Probes should be made with rprism, keeping in mind size limitations

import sys
import math
import numpy as np
import pdb

# Initializes probe as a rectangular prism outline of electrons with volume = 2*s3 * 2*s2 * 2*s1

def initProbe(x_c,y_c,xi_c,t0,s1,s2,xdensity,ydensity,xidensity,resolution):
    x_0, y_0, xi_0, z_0 = [],[],[],[]

    xistep = 2*s2/xidensity # Can also use resolution here
    ystep = 2*s1/ydensity
    xstep = xistep          # Purposefully setting xstep as equal to xistep for projection of x onto xi, unused here

# Calculate s3 
    s3 = xstep*(xdensity-1)/2.0
    print(f"{s3 = }")
    
# Define corners (front is first to enter field)
    ytop = y_c + s1
    ybot = y_c - s1
    xileft = xi_c - s2
    xiright = xi_c + s2
    xfront = x_c + s3  
    xback = x_c - s3

    print(f"Width of probe in x: 2*s3 = {2*s3}")
    print(f"Probe extends from x={xfront} to x={xback}")

# Start in front top right
    yn = ytop
    xin = xiright
    zn = xiright + t0
    xn = xfront

# Change to 2D Projection of x layering.
    xidensity_ = xidensity + xdensity - 1  # Allows enough particles in xi direction for x layering
    
# Place 2D particle sheet
    for i in range(0,ydensity):
        for j in range(0,xidensity_):
            x_0.append(xn)
            y_0.append(yn)
            xi_0.append(xin)
            z_0.append(xin + t0)
            xin -= xistep
        yn -= ystep
        xin = xiright

    return x_0, y_0, xi_0, z_0
