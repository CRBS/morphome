from __future__ import division

import math
import pyimod
import numpy as np

def transform_coords(model, x, y, z):
    """
    Transforms 3D coordinates from Amira to IMOD, relying upon the data entered
    in a model file's MINX chunk.
    """

    trans = [float(i) for i in [0, 0, 0]] 
    scale = [float(i) for i in [1, 1, 1]] 
    if model.minx_set:
        trans = model.minx_ctrans
        scale = model.minx_cscale 
    x = (x - trans[0]) / scale[0]
    y = (y - trans[1]) / scale[1]
    z = (z - trans[2]) / scale[2]
    return [x, y, z]

def parse_sav_csv(fname):
    """
    Parses a CSV file output from the 'Surface Area Volume' module of Amira.
    Extracts the number of triangles, surface area, and volume, and converts
    the surface area and volume from Angstroms to microns.
    """
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

def compute_sav_metrics(sa, vol):
    """
    Computes metrics from an object's surface area and volume.
    """
    pi = math.pi
    savRatio = sa / vol
    sphericity = (pi ** (1/3) * (6 * vol) ** (2/3)) / sa
    return savRatio, sphericity 

def is_outlier(points, thresh = 3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
