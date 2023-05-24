# For use with shape files of 'rprism_weighted_after' only


import sys
import math
import numpy as np
import pdb
import progressbar
import time

# Creates weights based on distribution and inputted masks (below)

def getWeights(beamx_c,beamy_c,beamxi_c,x_c,y_c,xi_c,s1,s2,xdensity,ydensity,xidensity,resolution,sigma_x,sigma_y,sigma_xi,noObj,t0,useWeights_x,useWeights_y,useWeights_xi,useMasks_x,useMasks_xi,useMasks_y):

# Recompute necessary parameters (as done in shape file)
    xidensity_ = xidensity + xdensity - 1  # Allows enough particles in xi direction for x layering
    xistep = 2*s2/xidensity # Can also use resolution here
    ystep = 2*s1/ydensity
    xstep = xistep          # Purposefully setting xstep as equal to xistep for projection of x onto xi

# Calculate s3 
    s3 = xstep*(xdensity-1)/2.0
    
# Define corners (xfront is first to enter field)
    ytop = y_c + s1
    ybot = y_c - s1
    xileft = xi_c - s2
    xiright = xi_c + s2
    xfront = x_c + s3  
    xback = x_c - s3
    
# Start in front top right
    yn = ytop
    xin = xiright
    xn = xfront
    zn = xiright + t0
    
    print("Creating weighting arrays...")
# Create individual coordinate arrays    
    x_0 = np.linspace(xfront,xback,xdensity)      # List of all possible inital x  positions of all particles going through simulation
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

    w_export1 = []
    w_export2 = w_y
    w_export3 = w_xi

# Loop through x layers to calculate weights with masks and add to 2D Y-Xi projection
    for i in progressbar.progressbar(range(0,len(x_0)), redirect_stout=False):
        start_time_weightcalc = time.time()

        if (useWeights_x):
            # Find weight value for this x-slice
            w_x = np.exp((-1.*(x_0[i]-beamx_c)**2)/(2*sigma_x**2)) # Gives a float

            # Check for x mask
            w_x = xMasks(useMasks_x,x_0[i],w_x)

            w_export1.append(w_x)

            w_y = w_y * w_x

        # Create final weighting list w to return
        # Maps 3d virtual particles in x-layer onto 2d Y-Xi projection appropriate location
        for k in range(0,len(xi_0)):
            # Multiply by xi weighting if in use
            w_virt = w_y * w_xi[k]
            
            # Loop through y layers and apply weighting in appropriate location
            for j in range(0,len(y_0)):
                w[xidensity_ * j + k + i] += w_virt[j]

        # Delete/Deallocate arrays for memory conservation
        w_x = None
        w_virt = None
        w_xy = None

    w_export1 = np.array(w_export1)
    w_export4 = w[int(xidensity_*(ydensity/2)):int(xidensity_*(ydensity/2)+xidensity_)]
    #np.savez("weight-exports.npz", w_exportx=w_export1, w_exporty=w_export2, w_exportxi=w_export3, w_exportvirt=w_export4)
    
    return w, w_export1, w_export2, w_export3, w_export4

def xiMasks(xi_0, w_xi):
    # Define masks in xi direction. Change if different mask is desired
    left_of_masks = [-26,-22]  # left most limit of each mask in order, on inital xi position
    right_of_masks = [-25,-20]  # right most limit of each mask in order, on initial xi position

    # Apply masks to w_xi
    for g in range(0,len(left_of_masks)):
        w_xi = np.where(np.logical_and(xi_0 > left_of_masks[g], xi_0 < right_of_masks[g]), 0, w_xi)

    return w_xi

def yMasks(y_0, w_y):
    # Define masks in y direction, 0 is 0 on the y-axis. Change if different mask is desired
    top_of_masks = [-0.3,0.4]  #upper limit of each mask in order, on inital y position
    bot_of_masks = [-0.4,0.2]  #lower limit of each mask in order, on inital y position

    # Apply masks to w_y
    for h in range(0,len(top_of_masks)):
        w_y = np.where(np.logical_and(y_0 > bot_of_masks[h], y_0 < top_of_masks[h]), 0, w_y)

    return w_y

def xMasks(useMasks_x,x_0_current,w_x):
    if (useMasks_x):
        # Define masks in x direction. Change if different mask is desired
        back_of_masks = []  # back limit of each mask in order, on inital x position
        front_of_masks = []  # right limit of each mask in order, on initial x position

        # Apply masks to w_x
        for m in range(0,len(back_of_masks)):
            if (np.logical_and(x_0_current > back_of_masks[m], x_0_current < front_of_masks[m])): # If in region of mask
                w_x = 0 # Set x weight to zero (0)

    return w_x