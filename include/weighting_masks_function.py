# For use with shape files of 'rectangle' only

import sys
import math
import numpy as np
import pdb
import progressbar
import time

# Creates weights based on distribution and inputted masks (below)
# 2D Version

def getWeights(beamx_c,beamy_c,beamxi_c,x_c,y_c,xi_c,s1,s2,xdensity,ydensity,xidensity,resolution,sigma_x,sigma_y,sigma_xi,noObj,t0,useWeights_x,useWeights_y,useWeights_xi,useMasks_x,useMasks_xi,useMasks_y):

# Recompute necessary parameters (as done in shape file)
    xistep = 2*s2/xidensity # Can also use resolution here
    ystep = 2*s1/ydensity
    
# Define corners (xfront is first to enter field)
    ytop = y_c + s1
    ybot = y_c - s1
    xileft = xi_c - s2
    xiright = xi_c + s2
    
# Start in front top right
    yn = ytop
    xin = xiright
    xn = x_c
    zn = xiright + t0
    
    print("Creating weighting arrays...")
# Create individual coordinate arrays    
    y_0 = np.linspace(ytop,ybot,ydensity)         # List of all possible inital y  positions of all particles going through simulation
    xi_0 = np.linspace(xiright,xileft,xidensity)  # List of all possible inital xi positions of all particles going through simulation

# Create empty weighting list
    w = []
    w = [0 for k in range(0,noObj)] # Creates weighting array of length noObj, with default value 0
    
# Masking and weighting for y and xi --------------------------------------------------------------------
    
    # Default w_y and w_xi weights
    w_y = np.full(y_0.shape, 1.0)
    w_xi = np.full(xi_0.shape, 1.0)

    if (useWeights_y): # If using gaussian weighting in y, apply to w_y
        w_y = np.exp((-1.*(y_0-beamy_c)**2)/(2*sigma_y**2)) # Calculate weights for each y slice
    
    if (useMasks_y): # If using Masks for y, apply them to y weighting array
        w_y = yMasks(y_0,w_y)
    
    if (useWeights_xi): # If using gaussian weighting in xi, apply to w_xi
        w_xi = np.exp((-1.*(xi_0-beamxi_c)**2)/(2*sigma_xi**2)) # Calculate weights for each xi slice

    if (useMasks_xi): # If using Masks for xi, apply them to xi weighting array
        w_xi = xiMasks(xi_0,w_xi)

    start_time_weightcalc = time.time()

    # Create final weighting list w to return
    for k in range(0,len(xi_0)):
        # Multiply by xi weighting if in use
        w_virt = w_y * w_xi[k]

        # Loop through y layers and apply weighting in appropriate location
        for j in range(0,len(y_0)):
            w[xidensity * j + k] += w_virt[j]

    
    return w, w[int(xidensity*39):int(xidensity*40-1)], w_y, w_xi


def xiMasks(xi_0, w_xi):
    # Define masks in xi direction. Change if different mask is desired
    left_of_masks = []  # left most limit of each mask in order, on inital xi position
    right_of_masks = []  # right most limit of each mask in order, on initial xi position

    # Apply masks to w_xi
    for g in range(0,len(left_of_masks)):
        w_xi = np.where(np.logical_and(xi_0 > left_of_masks[g], xi_0 < right_of_masks[g]), 0, w_xi)

    return w_xi

def yMasks(y_0, w_y):
    # Define masks in y direction, 0 is 0 on the y-axis. Change if different mask is desired
    top_of_masks = [0.5]#, 1.0]  #upper limit of each mask in order, on inital y position
    bot_of_masks = [0.4]#, 0.95]  #lower limit of each mask in order, on inital y position

    print(f"Masks extend from r={top_of_masks[0]} to r={bot_of_masks[0]}")
    # Apply masks to w_y
    for h in range(0,len(top_of_masks)):
        w_y = np.where(np.logical_and(y_0 > bot_of_masks[h], y_0 < top_of_masks[h]), 0, w_y)

    return w_y