# This file retrieves the fields from OSIRIS Quasi3D data files stored within the data/ folder.
# This specific version is an adaptation from previous programs that read data from .h5 files produced by UCLA. This program will read in .h5 files created by Arohi Jain at Stony Brook.
# Note that these files do NOT consider spherical harmonics of the laser; it is just EM fields of wake and fields of Laser. 
# Functions that MUST be updated for any new simulation (i.e. Called in either main.py or eProbe.py) are designated with three asterisks ***
# All other functions are used for either reading out data or plotting results

import sys
import h5py as h5
import numpy as np
import math
import pdb

# Coordinate System
# z   - Direction of laser propagation (longitudinal)
# xi  - Position along z relative to wavefront
# r   - Cylindrical coordinate around z
# phi - Cylindrical coordinate around z, define phi = 0 along x
# x   - Direction of transverse probe
# y   - Direction perpendicular to transverse probe

# Modes
# mode = 0 refers to LWF effects only
# mode = 1 refers to laser effects only
# mode = any other integer uses LWF + laser effects

# Definition of Constants
M_E = 9.109e-31                       # Electron rest mass in kg
EC = 1.60217662e-19                   # Electron charge in C
EP_0 = 8.854187817e-12                # Vacuum permittivity in C/(V m)
C = 299892458                         # Speed of light in vacuum in m/s

Quasi_ID = '10e14' #field data for 10e14 plasma density

# order of entries in dataset: 
#'By_L', 'By_wake', 'Ex_L', 'Ex_wake', 'Ez_L', 'Ez_wake', 'r', 'xi'

def getFile():
    
    fpath = 'data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5'
    #fpath = 'data/OSIRIS/Quasi3D/b3_cyl_m-1-im-'+ Quasi_ID + '.h5'

    return fpath
    
    
def getField(fpath): 
    
    f = h5.File(fpath,"r")
    #print(f"{f}")
    datasetNames = [n for n in f.keys()]
    #fields = datasetNames[0:6]
    #print(fields)
    Field_dat = []
    #print(type(f[fields]), len(f[fields]))
    #print(f[fields[:]][:])
    for i in range(len(datasetNames)):
        Field_dat.append(f[datasetNames[i]][:].astype(float))
    
    #print(len(Field_dat))
    
    return Field_dat

#getField('data/OSIRIS/Quasi3D/b3_cyl_m-1-im-'+ Quasi_ID + '.h5')
"""
def getTime(): # ***
    f = h5.File('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5',"r")
    t0 = f.attrs['TIME']
    t0 = t0[0]
    return t0
"""
def getPlasDensity(): 
    if (Quasi_ID == '10e14'):
        return 1e15
    else:
        return 10e15

def getPropagationSpeed(): # Define the group velocity of the laser
    if (Quasi_ID == '10e14'):
        return 1.000 #THIS IS INCORRECT, JUST FOR TESTING
    #elif (Quasi_ID == '000067'):
        return 0.9958959
    else:
        return 1

def getPlasFreq(): 
    N_0 = getPlasDensity()
    return math.sqrt(EC**2 * N_0 / (M_E * EP_0))

def getFieldNormal(): 
    
    W_p = getPlasFreq()
    Enorm = M_E*C*W_p/EC 
    
    return Enorm

def axes(): 
# Retrieve axes boundaries under staggered mesh
    f = h5.File('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5',"r")
    datasetNames = [n for n in f.keys()] # Three Datasets: AXIS, SIMULATION, Field data
    #field = datasetNames[-2:0]
    #Field_dat = f[field][:].astype(float)
    xi_data = f['xi'] # xi-points, xi[0] is xi_max
    r_data = f['r'] # radial points, r[-1] is r_max
    dz = 2.082e-3 # Width of each cell
    dr = 8.921e-4
    
# Account for specific definitions of bottom-left values for each field component
    xi_bounds_1 = [xi_data[-1], xi_data[0]] # Used for E1, B1
    #xi_bounds_2 = [a1_bounds[0] - dz/2, a1_bounds[1] + dz/2] # Used for E2, E3, B2, B3
    r_bounds_1 = [r_data[0], r_data[-1]] # Used for E2, B2
    #r_bounds_2 = [a2_bounds[0] - dr/2, a2_bounds[1] + dr] # Used for E1, E3, B1, B3

# Field Shape is (433, 25231), where data is written as E(r,z)
    xiaxis_1 = np.linspace(xi_bounds_1[0], xi_bounds_1[1], len(xi_data))
    raxis_1 = np.linspace(r_bounds_1[0], r_bounds_1[1], len(r_data))

    return xi_data, r_data, xiaxis_1, raxis_1

rpoints, xipoints, xiaxis_1, raxis_1 = axes() # Evenly spaced axes data
"""
def getBoundCond(): # ***
# Define when the electron leaves the plasma cell
    f = h5.File('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5',"r")
    datasetNames = [n for n in f.keys()] # Three Datasets: AXIS, SIMULATION, Field data
    field = datasetNames[-1]
    Field_dat = f[field][:].astype(float)
    a1_bounds = f['AXIS']['AXIS1'] # zmin and zmax
    a2_bounds = f['AXIS']['AXIS2'] # rmin and rmax - dr/2
    dr = 6.59e-3
    t0 = getTime() # Call getTime func to get sim time
    return [a1_bounds[0] - t0, a1_bounds[1] - t0, a2_bounds[1] + dr/2] # ximin, ximax, rmax

# Return cylindrical Electric field components
# E1 - z
# E2 - r
# E3 - phi
"""

def getWakeEMFields(): 
    
    wake_By = getField('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5')[1]
    wake_Ex = getField('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5')[3]
    wake_Ez = getField('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5')[5]
    
    return wake_By, wake_Ex, wake_Ez

def getLaserEMFields(): 
    
    laser_By = getField('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5')[0]
    laser_Ex = getField('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5')[2]
    laser_Ez = getField('data/OSIRIS/Quasi3D/fields_'+ Quasi_ID + '.h5')[4]
    
    return laser_By, laser_Ex, laser_Ez


"""
def getchargeElectrons_M0():
    return getField('data/OSIRIS/Quasi3D/charge_cyl_m-electrons-0-re-'+ Quasi_ID + '.h5')

def getchargeElectrons_M1_Re():
    return getField('data/OSIRIS/Quasi3D/charge_cyl_m-electrons-1-re-'+ Quasi_ID + '.h5')

def getchargeElectrons_M1_Im():
    return getField('data/OSIRIS/Quasi3D/charge_cyl_m-electrons-1-im-'+ Quasi_ID + '.h5')

def getchargeIons_M0():
    return getField('data/OSIRIS/Quasi3D/charge_cyl_m-ions-0-re-'+ Quasi_ID + '.h5')

def getchargeIons_M1_Re():
    return getField('data/OSIRIS/Quasi3D/charge_cyl_m-ions-1-re-'+ Quasi_ID + '.h5')

def getchargeIons_M1_Im():
    return getField('data/OSIRIS/Quasi3D/charge_cyl_m-ions-1-im-'+ Quasi_ID + '.h5')
"""

wakeBy, wakeEx, wakeEz = getWakeEMFields()
laserBy, laserEx, laserEz = getLaserEMFields()

"""
if Quasi_ID=="000067":
    chargeElectrons_M0 = getchargeElectrons_M0()
    chargeElectrons_M1_Re = getchargeElectrons_M1_Re()
    chargeElectrons_M1_Im = getchargeElectrons_M1_Im()
    chargeIons_M0 = getchargeIons_M0()
    chargeIons_M1_Re = getchargeIons_M1_Re()
    chargeIons_M1_Im = getchargeIons_M1_Im()


def getPhi(x,y):
    return math.atan2(y,x) # From -pi to pi

def find_nearest_index(array,value):
    idx = np.searchsorted(array, value, side="right")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def EField(axis,x,y,xi,r,vx=-1,vy=-1,vz=-1,vr=-1,vphi=-1,mode=-1): # ***
# axis = 1 refers to z-axis field
# axis = 2 refers to x-axis field
# axis = 3 refers to y-axis field
# mode = 0 refers to LWF effects only
# mode = 1 refers to laser effects only
# mode = any other integer uses LWF + laser effects
    phi = getPhi(x,y)
    cos = math.cos(phi)
    sin = math.sin(phi)
    xiDex1 = find_nearest_index(xiaxis_1, xi)
    xiDex2 = find_nearest_index(xiaxis_2, xi)
    rDex1 = find_nearest_index(raxis_1, r)
    rDex2 = find_nearest_index(raxis_2, r)
    # Return expanded EFields
    if (mode == 0):
        if (axis == 1):
            return E1_M0[rDex2, xiDex1]
        elif (axis == 2):
            return E2_M0[rDex1, xiDex2]*cos - E3_M0[rDex2, xiDex2]*sin
        elif (axis == 3):
            return E3_M0[rDex2, xiDex2]*cos + E2_M0[rDex1, xiDex2]*sin
    elif (mode == 1):
        if (axis == 1):
            return E1_M1_Re[rDex2, xiDex1]*cos + E1_M1_Im[rDex2, xiDex1]*sin
        elif (axis == 2):
            return E2_M1_Re[rDex1, xiDex2]*cos**2 - E3_M1_Re[rDex2, xiDex2]*cos*sin + E2_M1_Im[rDex1, xiDex2]*cos*sin - E3_M1_Im[rDex2, xiDex2]*sin**2
        elif (axis == 3):
            return E3_M1_Re[rDex2, xiDex2]*cos**2 + E2_M1_Re[rDex1, xiDex2]*cos*sin + E3_M1_Im[rDex2, xiDex2]*cos*sin + E2_M1_Im[rDex1, xiDex2]*sin**2
    else:
        if (axis == 1):
            return E1_M0[rDex2, xiDex1] + E1_M1_Re[rDex2, xiDex1]*cos + E1_M1_Im[rDex2, xiDex1]*sin
        elif (axis == 2):
            return E2_M0[rDex1, xiDex2]*cos - E3_M0[rDex2, xiDex2]*sin + E2_M1_Re[rDex1, xiDex2]*cos**2 - E3_M1_Re[rDex2, xiDex2]*cos*sin + E2_M1_Im[rDex1, xiDex2]*cos*sin - E3_M1_Im[rDex2, xiDex2]*sin**2
        elif (axis == 3):
            return E3_M0[rDex2, xiDex2]*cos + E2_M0[rDex1, xiDex2]*sin + E3_M1_Re[rDex2, xiDex2]*cos**2 + E2_M1_Re[rDex1, xiDex2]*cos*sin + E3_M1_Im[rDex2, xiDex2]*cos*sin + E2_M1_Im[rDex1, xiDex2]*sin**2

def BForce(axis,x,y,xi,r,vx=-1,vy=-1,vz=-1,vr=-1,vphi=-1,mode=-1): # ***
# axis = 1 refers to z-axis field
# axis = 2 refers to x-axis field
# axis = 3 refers to y-axis field
    phi = getPhi(x,y)
    cos = math.cos(phi)
    sin = math.sin(phi)
    xiDex1 = find_nearest_index(xiaxis_1, xi)
    xiDex2 = find_nearest_index(xiaxis_2, xi)
    rDex1 = find_nearest_index(raxis_1, r)
    rDex2 = find_nearest_index(raxis_2, r)
    # Calculate expanded BFields
    if (mode == 0):
        Bz = B1_M0[rDex2, xiDex1]
        Bx = B2_M0[rDex1, xiDex2]*cos - B3_M0[rDex2, xiDex2]*sin
        By = B3_M0[rDex2, xiDex2]*cos + B2_M0[rDex1, xiDex2]*sin
    elif (mode == 1):
        Bz = B1_M1_Re[rDex2, xiDex1]*cos + B1_M1_Im[rDex2, xiDex1]*sin
        Bx = B2_M1_Re[rDex1, xiDex2]*cos**2 - B3_M1_Re[rDex2, xiDex2]*cos*sin + B2_M1_Im[rDex1, xiDex2]*cos*sin - B3_M1_Im[rDex2, xiDex2]*sin**2
        By = B3_M1_Re[rDex2, xiDex2]*cos**2 + B2_M1_Re[rDex1, xiDex2]*cos*sin + B3_M1_Im[rDex2, xiDex2]*cos*sin + B2_M1_Im[rDex1, xiDex2]*sin**2
    else:
        Bz = B1_M0[rDex2, xiDex1] + B1_M1_Re[rDex2, xiDex1]*cos + B1_M1_Im[rDex2, xiDex1]*sin
        Bx = B2_M0[rDex1, xiDex2]*cos - B3_M0[rDex2, xiDex2]*sin + B2_M1_Re[rDex1, xiDex2]*cos**2 - B3_M1_Re[rDex2, xiDex2]*cos*sin + B2_M1_Im[rDex1, xiDex2]*cos*sin - B3_M1_Im[rDex2, xiDex2]*sin**2
        By = B3_M0[rDex2, xiDex2]*cos + B2_M0[rDex1, xiDex2]*sin + B3_M1_Re[rDex2, xiDex2]*cos**2 + B2_M1_Re[rDex1, xiDex2]*cos*sin + B3_M1_Im[rDex2, xiDex2]*cos*sin + B2_M1_Im[rDex1, xiDex2]*sin**2
    # Cross-product velocities with BFields and return the BForce
    if axis == 1:
        return vx * By - vy * Bx
    elif axis == 2:
        return vy * Bz - vz * By
    elif axis == 3:
        return -1.0 * (vx * Bz - vz * Bx)

def BField(axis,x,y,xi,r,vx=-1,vy=-1,vz=-1,vr=-1,vphi=-1,mode=-1):
# Return BField
    phi = getPhi(x,y)
    cos = math.cos(phi)
    sin = math.sin(phi)
    xiDex1 = find_nearest_index(xiaxis_1, xi)
    xiDex2 = find_nearest_index(xiaxis_2, xi)
    rDex1 = find_nearest_index(raxis_1, r)
    rDex2 = find_nearest_index(raxis_2, r)
# Calculate expanded BFields
    if (mode == 0):
        if (axis == 1):
            return B1_M0[rDex2, xiDex1]
        elif (axis == 2):
            return B2_M0[rDex1, xiDex2]*cos - B3_M0[rDex2, xiDex2]*sin
        elif (axis == 3):
            return B3_M0[rDex2, xiDex2]*cos + B2_M0[rDex1, xiDex2]*sin
    elif (mode == 1):
        if (axis == 1):
            return B1_M1_Re[rDex2, xiDex1]*cos + B1_M1_Im[rDex2, xiDex1]*sin
        elif (axis == 2):
            return B2_M1_Re[rDex1, xiDex2]*cos**2 - B3_M1_Re[rDex2, xiDex2]*cos*sin + B2_M1_Im[rDex1, xiDex2]*cos*sin - B3_M1_Im[rDex2, xiDex2]*sin**2
        elif (axis == 3):
            return B3_M1_Re[rDex2, xiDex2]*cos**2 + B2_M1_Re[rDex1, xiDex2]*cos*sin + B3_M1_Im[rDex2, xiDex2]*cos*sin + B2_M1_Im[rDex1, xiDex2]*sin**2
    else:
        if (axis == 1):
            return B1_M0[rDex2, xiDex1] + B1_M1_Re[rDex2, xiDex1]*cos + B1_M1_Im[rDex2, xiDex1]*sin
        elif (axis == 2):
            return B2_M0[rDex1, xiDex2]*cos - B3_M0[rDex2, xiDex2]*sin + B2_M1_Re[rDex1, xiDex2]*cos**2 - B3_M1_Re[rDex2, xiDex2]*cos*sin + B2_M1_Im[rDex1, xiDex2]*cos*sin - B3_M1_Im[rDex2, xiDex2]*sin**2
        elif (axis == 3):
            return B3_M0[rDex2, xiDex2]*cos + B2_M0[rDex1, xiDex2]*sin + B3_M1_Re[rDex2, xiDex2]*cos**2 + B2_M1_Re[rDex1, xiDex2]*cos*sin + B3_M1_Im[rDex2, xiDex2]*cos*sin + B2_M1_Im[rDex1, xiDex2]*sin**2


def chargeElectrons(axis,x,y,xi,r,vx=-1,vy=-1,vz=-1,vr=-1,vphi=-1,mode=-1): # ***
# axis = 1 refers to z-axis field
# axis = 2 refers to x-axis field
# axis = 3 refers to y-axis field
# mode = 0 refers to LWF effects only
# mode = 1 refers to laser effects only
# mode = any other integer uses LWF + laser effects
    phi = getPhi(x,y)
    cos = math.cos(phi)
    sin = math.sin(phi)
    xiDex1 = find_nearest_index(xiaxis_1, xi)
    xiDex2 = find_nearest_index(xiaxis_2, xi)
    rDex1 = find_nearest_index(raxis_1, r)
    rDex2 = find_nearest_index(raxis_2, r)
    # Return chargeElectrons
    if (mode == 0):
        if (axis == 1):
            return chargeElectrons_M0[rDex2, xiDex1]
    elif (mode == 1):
        if (axis == 1):
            return chargeElectrons_M1_Re[rDex2, xiDex1]*cos + chargeElectrons_M1_Im[rDex2, xiDex1]*sin
    else:
        if (axis == 1):
            return chargeElectrons_M0[rDex2, xiDex1] + chargeElectrons_M1_Re[rDex2, xiDex1]*cos + chargeElectrons_M1_Im[rDex2, xiDex1]*sin


def chargeIons(axis,x,y,xi,r,vx=-1,vy=-1,vz=-1,vr=-1,vphi=-1,mode=-1): # ***
# axis = 1 refers to z-axis field
# axis = 2 refers to x-axis field
# axis = 3 refers to y-axis field
# mode = 0 refers to LWF effects only
# mode = 1 refers to laser effects only
# mode = any other integer uses LWF + laser effects
    phi = getPhi(x,y)
    cos = math.cos(phi)
    sin = math.sin(phi)
    xiDex1 = find_nearest_index(xiaxis_1, xi)
    xiDex2 = find_nearest_index(xiaxis_2, xi)
    rDex1 = find_nearest_index(raxis_1, r)
    rDex2 = find_nearest_index(raxis_2, r)
    # Return chargeIons
    if (mode == 0):
        if (axis == 1):
            return chargeIons_M0[rDex2, xiDex1]
    elif (mode == 1):
        if (axis == 1):
            return chargeIons_M1_Re[rDex2, xiDex1]*cos + chargeIons_M1_Im[rDex2, xiDex1]*sin
    else:
        if (axis == 1):
            return chargeIons_M0[rDex2, xiDex1] + chargeIons_M1_Re[rDex2, xiDex1]*cos + chargeIons_M1_Im[rDex2, xiDex1]*sin
"""    