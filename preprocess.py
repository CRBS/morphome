import pyimod

def nucleus(fnamein, fnameout):
    mod = pyimod.ImodModel(fnamein)

    # Remove small contours (i.e. contours having <= 2 points. These are likely
    # erroneous, and are prone to creating spikes in the surface in Amira that
    # result in NaNs at vertices when computing the Shape Index.
    mod.removeSmallContours()
    mod.Objects[0].sortContours()

    # Perform initial smoothing of the surface using IMOD. This is primarily to
    # eliminate the effects of bad alignment and warping. Done in 3 steps:
    #    1. Remove pre-existing meshs
    #    2. Smooth the contours in 3D
    #    3. Remesh using a low resolution mesh 
    mod = pyimod.ImodCmd(mod, 'imodmesh -e')
    mod = pyimod.ImodCmd(mod, 'smoothsurf -nz 15 -dist 25')
    mod = pyimod.ImodCmd(mod, 'imodmesh -CT -l')

    pyimod.ImodExport(mod, fnameout, object = 1)
