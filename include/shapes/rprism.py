import sys
import math
import numpy as np
import pdb

# Initializes probe as a rectangular prism outline of electrons with volume = 2*s3 * 2*s2 * 2*s1

def initProbe(x_c,y_c,xi_c,t0,s1,s2,xdensity,ydensity,xidensity,resolution):
    x_0, y_0, xi_0, z_0 = [],[],[],[]

    xistep = 2*s2/xidensity # Can also use resolution here
    ystep = 2*s1/ydensity
    xstep = xistep

    s3 = xstep*(xdensity-1)/2.0


# Define corners (front is first to enter field)
    ytop = y_c + s1
    ybot = y_c - s1
    xileft = xi_c - s2
    xiright = xi_c + s2
    xfront = x_c + s3  
    xback = x_c - s3

# Start in front top left
    yn = ytop
    xin = xileft
    zn = xileft + t0
    xn = xfront

    for i in range(0,xdensity):
        for j in range(0,ydensity):
            for k in range(0,xidensity):
                x_0.append(xn)
                y_0.append(yn)
                xi_0.append(xin)
                z_0.append(xin + t0)
                xin += xistep
            yn -= ystep
            xin = xileft
        print('x Layer',i+1,'of',xdensity,'complete') #Just for debugging // can remove
        xn -= xstep
        yn = ytop

    return x_0, y_0, xi_0, z_0
