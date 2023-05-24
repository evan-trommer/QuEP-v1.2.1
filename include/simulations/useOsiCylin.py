# This file retrieves the fields from cylindrically symmetric OSIRIS data files stored within the data/ folder.

import h5py as h5
import numpy as np
import math
import pdb

# Coordinate System
# z   - Direction of laser propagation (longitudinal)
# xi  - Position along z relative to wavefront (unnecessary?)
# r   - Cylindrical coordinate around z
# phi - Cylindrical coordinate around z, define phi = 0 along x
# x   - Direction of transverse probe
# y   - Direction perpendicular to transverse probe

def getField(fname):
  f = h5.File(fname,"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  return Field_dat

def axes():
  f = h5.File('data/OSIRIS/CylinSymm/EField_r.h5',"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  a1_bounds = f['AXIS']['AXIS1']
  a2_bounds = f['AXIS']['AXIS2']

  t0 = getTime() #time at which field data was simulated, constant for all fields
  xi_dat = np.linspace(a1_bounds[0] - t0,a1_bounds[1] - t0,len(Field_dat[0]))
  r_dat = np.linspace(a2_bounds[0],a2_bounds[1],len(Field_dat))

  return r_dat, xi_dat

def getTime():
    return 858.95 # Time at which field data was simulated, constant for all fields

def getPlasDensity():
    return 1e23
    
def transE():
  return getField('data/OSIRIS/CylinSymm/EField_r.h5')

def longE():
  return getField('data/OSIRIS/CylinSymm/EField_z.h5')

def phiB():
  return getField('data/OSIRIS/CylinSymm/BField_phi.h5')

r_sim, xi_sim = axes()
Er_sim = transE()
Ez_sim = longE()
Bphi_sim = phiB()

def getPhi(x,y):
    return math.atan2(y,x) # From -pi to pi

def find_nearest_index(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def EField(axis,x,y,z,r=-1,vx=-1,vy=-1,vz=-1,vr=-1,vphi=-1):
# axis = 2 refers to x-axis field
# axis = 3 refers to y-axis field
# axis = 1 refers to z-axis field
# axis = 4 refers to r-axis field
    phi = getPhi(x, y)
    zDex = find_nearest_index(xi_sim, z)
    rDex = find_nearest_index(r_sim, r)
    if axis == 2:
        return Er_sim[rDex, zDex] * math.cos(phi)
    elif axis == 3:
        return Er_sim[rDex, zDex] * math.sin(phi)
    elif axis == 1:
        return Ez_sim[rDex, zDex]
    elif axis == 4:
        return -1.0*Er_sim[rDex, zDex]

def BForce(axis,x,y,z,r=-1,vx=-1,vy=-1,vz=-1,vr=-1,vphi=-1):
# axis = 2 refers to x-axis field
# axis = 3 refers to y-axis field
# axis = 1 refers to z-axis field
    phi = getPhi(x, y)
    zDex = find_nearest_index(xi_sim, z)
    rDex = find_nearest_index(r_sim, r)
    if axis == 2:
        return -1.0 * vz * Bphi_sim[rDex, zDex] * math.cos(phi)
    elif axis == 3:
        return -1.0 * vz * Bphi_sim[rDex, zDex] * math.sin(phi)
    elif axis == 1:
        return vx * Bphi_sim[rDex, zDex] * math.cos(phi) + vy * Bphi_sim[rDex, zDex] * math.sin(phi)

def getBoundCond():
    t0 = getTime()
    return [858 - t0, 868 - t0, 6]
