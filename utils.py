import pyimod

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
