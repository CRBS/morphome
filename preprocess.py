import pyimod

def general_prep(fnamein):
    # If a file is supplied, load the file. If not, proceed with the ImodModel
    # object.
    typein = type(fnamein).__name__
    if typein == 'str':
        mod = pyimod.ImodModel(fnamein)
    elif typein == 'ImodModel':
        mod = fnamein

    # Remove small contours (i.e. contours having <= 2 points. These are likely
    # erroneous, and are prone to creating spikes in the surface in Amira that
    # can result in NaNs at vertices. 
    mod.removeSmallContours()
    mod.Objects[0].sortContours()

    # Remove pre-existing meshes
    mod = pyimod.ImodCmd(mod, 'imodmesh -e')

    return mod   

def mitochondrion(fnamein, fnameout)
    mod = general_prep(fnamein)
    # Perform initial smoothing of the surface using IMOD. Done in 2 steps:
    #    1. Smooth the contours in 3D
    #    2. Remesh. Do not use a low resolution mesh, as this tends to contract
    #       mitochondria.
    mod = pyimod.ImodCmd(mod, 'smoothsurf -nz 15 -dist 25')
    mod = pyimod.ImodCmd(mod, 'imodmesh -CTs -P 4')
    pyimod.ImodExport(mod, fnameout, object = 1)

def nucleus(fnamein, fnameout):
    mod = general_prep(fnamein)
    # Perform initial smoothing of the surface using IMOD. This is primarily to
    # eliminate the effects of bad alignment and warping. Done in 2 steps:
    #    1. Smooth the contours in 3D
    #    2. Remesh using a low resolution mesh 
    mod = pyimod.ImodCmd(mod, 'smoothsurf -nz 15 -dist 25')
    mod = pyimod.ImodCmd(mod, 'imodmesh -CT -l')
    pyimod.ImodExport(mod, fnameout, object = 1)
