#!/usr/bin/env python
from __future__ import division 
import numpy as np
import re

"""
Takes a vertex index, or seed, as input, then finds all triangles that
contain this seed vertex. The vertex indices of all these neighboring
triangles are then returned. 
"""
def get_adjacent_vertices(seed, indices):
    # Find all triangles that contain the seed vertex
    tris_with_seed = np.where(indices == seed)[0]

    # Find all vertices of these triangles. Keep only the unique ones
    idx = []
    for i in tris_with_seed:
        idx = np.hstack([idx, indices[i]])
    idx = np.unique(idx).tolist()
    idx = [int(x) for x in idx]
    idx.remove(seed)
    return idx

"""

"""
def find_nuclear_folds(model, vertices, indices, shapeIndex, grad):
    coordList = []
    nTri = indices.shape[0]
    nVert = vertices.shape[0]
    idx_keep = []

    # Find all vertices with Shape Index > 0
    vert_idx_pos = np.where(shapeIndex > 0)[0].tolist()
    vert_idx_pos = [int(x) + 1 for x in vert_idx_pos]

    nCluster = 0
    while vert_idx_pos:
        # Set the seed vertex as the first entry
        vert_seed = vert_idx_pos[0]
        vert_check = [vert_seed]
        vert_keep = []

        while vert_check:
            # Find all vertices in triangles adjacent to the seed vertex
            vert_adj = get_adjacent_vertices(vert_check[0], indices)
            vert_keep.append(vert_check[0])

            # If none of the adjacent vertices have a positive Shape Index,
            # then the seed vertex is an island. Break the loop and remove it.
            # If there is at least one adjacent vertex with a positive Shape
            # Index, add the seed vertex to the keep list and remove it from
            # the seed position of the check list.
            if not sum([x in vert_idx_pos for x in vert_adj]):
                break
            else:
                vert_check = vert_check[1:]

            # Append all adjacent vertices with positive Shape Indices to the
            # check list.
            for i in vert_adj:
                if ((i in vert_idx_pos) and (i not in vert_check) and
                    (i not in vert_keep)):
                    vert_check.append(i)

        for i in vert_keep:
            vert_idx_pos.remove(i)

        # Process the relevantly sized clusters
        if len(vert_keep) > 20:
            nCluster+=1
            nVertCluster = len(vert_keep)

            # Get a list of Shape Index values from the vertices
            K = 0
            sCluster = np.zeros([nVertCluster, 1]) 
            for i in vert_keep:
                sCluster[K] = shapeIndex[i - 1]
                K+=1

            # Get the triangles from vertices
            tri_keep = []
            for i in vert_keep:
                tris_with_vert = np.where(indices == i)[0]
                for j in tris_with_vert:
                    idx = indices[j]
                    if (sum([x in vert_keep for x in idx]) == 3 and
                        j not in tri_keep):
                        tri_keep.append(j)
            nTriCluster = len(tri_keep)

            # Get the gradients from the triangles. Take the magnitude of the
            # Shape Index gradient at each triangle face.
            magMinCluster = 1
            for i in tri_keep:
                gradi = grad[i,:]
                magClusteri = np.linalg.norm(gradi)
                if magClusteri < magMinCluster:
                    magMinCluster = magClusteri
                    imc = indices[i,:]

            # Calculate the coordinates of the triangle centroid
            coordx = (vertices[imc[0]-1,0] + vertices[imc[1]-1,0] +
                vertices[imc[2]-1,0]) / 3
            coordy = (vertices[imc[0]-1,1] + vertices[imc[1]-1,1] +
                vertices[imc[2]-1,1]) / 3
            coordz = (vertices[imc[0]-1,2] + vertices[imc[1]-1,2] +
                vertices[imc[2]-1,2]) / 3   
           
            # Scale back to IMOD coordinate system
            trans = [float(x) for x in [0, 0, 0]]
            scale = [float(x) for x in [1, 1, 1]]
            if model.minx_set:
                trans = model.minx_ctrans
                scale = model.minx_cscale 
            coordx = (coordx + trans[0]) / scale[0]
            coordy = (coordy + trans[1]) / scale[1]
            coordz = (coordz + trans[2]) / scale[2]

            # Print metrics
            print "Cluster {0}".format(nCluster)
            print "==========="
            print "# Vertices : {0}".format(nVertCluster)
            print "# Triangles: {0}".format(nTriCluster)
            print "S_min  : {0}".format(np.amin(sCluster))
            print "S_max  : {0}".format(np.amax(sCluster))
            print "S_mean : {0}".format(np.mean(sCluster))
            print "S_std  : {0}".format(np.std(sCluster))
            print "Min. gradient magnitude: {0}".format(magMinCluster)
            print "Coords: {0}, {1}, {2}".format(coordx, coordy, coordz)
            print ""
            
            coordList.append([coordx, coordy, coordz])
    return coordList
