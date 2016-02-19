#!/usr/bin/env python
from __future__ import division 
import numpy as np
import re

"""
Reads an Amira Mesh ASCII file of scalar field values to a Numpy array
"""
def scalar_field(fname):
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

"""
Reads an Amira Mesh ASCII file of vector field directions to a Numpy array
"""
def vector_field(fname):
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

"""
Read an HxSurface ASCII file, which describes (1) the mesh vertices and
(2) the vertex indices corresponding to all mesh triangles.
"""
def surface(fname):
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
