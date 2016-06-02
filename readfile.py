#!/usr/bin/env python
from __future__ import division 
import numpy as np
import re
import math

def scalar_field(fname):
    """
    Reads an Amira Mesh ASCII file of scalar field values to a Numpy array
    """
    C = 0 
    K = 0 
    with open(fname, 'r') as fid:
        for line in fid:
            if re.match('^nNodes', line) or re.match('^nTriangles', line):
                nVert = int(line.split()[1])
                sField = np.zeros([nVert, 1]) 
            if re.match('^@1', line):
                C = 1 
                continue
            if C and line.strip():
                sField[K] = float(line)
                K+=1
    return sField

def vector_field(fname):
    """
    Reads an Amira Mesh ASCII file of vector field directions to a Numpy array
    """
    C = 0
    K = 0
    with open(fname, 'r') as fid:
        for line in fid:
            if re.match('^nTriangles', line):
                nTri = int(line.split()[1])
                vField = np.zeros([nTri, 3])
            if re.match('^@1', line):
                C = 1
                continue
            if C and line.strip():
                vector = [np.absolute(float(x)) for x in line.split()]
                vField[K,] = vector
                K+=1
    return vField

def surface(fname):
    """
    Read an HxSurface ASCII file, which describes (1) the mesh vertices and
    (2) the vertex indices corresponding to all mesh triangles.
    """

    C = 0
    with open(fname, 'r') as fid:
        for line in fid: 
            if re.match('^Vertices', line):
                nVert = int(line.split()[1])
                C = 1
                Kvert = 0
                continue
            if re.match('^Triangles', line):
                nInd = int(line.split()[1])
                C = 2
                Kind = 0
                continue
            if C == 1 and Kvert < nVert:
                coords = [float(x) for x in line.split()]
                try:
                    vertices
                except NameError:
                    vertices = np.array(coords)
                else:
                    vertices = np.vstack([vertices, coords])   
                Kvert+=1 
            elif C == 1 and Kvert >= nVert:
                C = 0
                continue
            if C == 2 and Kind < nInd:
                inds = [int(x) for x in line.split()]
                try:
                    indices
                except NameError:
                    indices = np.array(inds)
                else:
                    indices = np.vstack([indices, inds])
                Kind+=1
            elif C == 2 and Kind >= nInd:
                C = 0
                continue
        return vertices, indices

def label_csv(fname):
    """
    Reads data from a CSV file output by Amira's 'Label Analysis' module. Data
    are returned as a list.
    """
    fid = open(fname, 'r')
    for i, line in enumerate(fid):
       if i == 2:
           metrics = line
           break
    metricList = line.split(',')
    metricList = metricList[0:-2]
    metricList = [float(x) for x in metricList]
    return metricList

def skel_length(fname):
    """
    Reads data from a CSV file output by Amira's spatial graph module.
    """
    with open(fname, 'r') as fid:
        for line in fid:
            if 'Graph Summary' in line:
                fid.next()
                vals = fid.next().split(',')
                nSegments = int(vals[1])
                lengthTot = float(vals[5])
                segments = np.zeros([nSegments, 4])
            elif 'Segment Statistics' in line:
                fid.next()
                for i in range(nSegments):
                    vals = fid.next().split(',')
                    segments[i,0] = float(vals[1])
                    segments[i,1] = float(vals[2])
                    segments[i,2] = float(vals[4])
                    segments[i,3] = float(vals[5])
            elif 'Node Statistics' in line:
                fid.next()
                nodes = [int(x) for x in fid.next().split(',')]
    print nSegments, lengthTot
    print segments
    print nodes

def sav(fname):
    """
    Reads data from a CSV file output by Amira's Surface Area module. Returns
    volume and surface area.
    """
    with open(fname, 'r') as fid:
        for line in fid:
            if 'Exterior' in line:
                line = line.split(',')
                sa = line 
                sa = math.fabs(int(line[2]))
                volume = math.fabs(int(line[3]))
                break
    sa = sa / (10000 ** 2)
    volume = volume / (10000 ** 3)
    return sa, volume
