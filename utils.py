from __future__ import division

import math
import pyimod

"""
Transforms 3D coordinates from Amira to IMOD, relying upon the data entered
in a model file's MINX chunk.
"""
def transform_coords(model, x, y, z):
    trans = [float(i) for i in [0, 0, 0]] 
    scale = [float(i) for i in [1, 1, 1]] 
    if model.minx_set:
        trans = model.minx_ctrans
        scale = model.minx_cscale 
    x = (x - trans[0]) / scale[0]
    y = (y - trans[1]) / scale[1]
    z = (z - trans[2]) / scale[2]
    return [x, y, z]

"""
Parses a CSV file output from the 'Surface Area Volume' module of Amira.
Extracts the number of triangles, surface area, and volume, and converts the
surface area and volume from Angstroms to microns.
"""
def parse_sav_csv(fname):
    fid = open(fname, 'r')
    for i, line in enumerate(fid):
       if i == 2:
           line = line.split(',')
           nTri = int(line[1])
           sa = abs(int(line[2]))
           vol = abs(int(line[3]))
           break
    fid.close()
    sa = sa / (10000 ** 2)
    vol = vol / (10000 ** 3) 
    return nTri, sa, vol

"""
Computes metrics from an object's surface area and volume
"""
def compute_sav_metrics(sa, vol):
    pi = math.pi
    savRatio = sa / vol
    sphericity = (pi ** (1/3) * (6 * vol) ** (2/3)) / sa
    return savRatio, sphericity 
