#!/usr/bin/env python

import sys
import os
import re
import shutil
import subprocess
import datetime
import pyimod
import morphome
import numpy as np
from optparse import OptionParser

def parse_args():
    global p
    p = OptionParser(usage = "%prog [options] file.mod /path/for/output")
    (opts, args) = p.parse_args()
    fileModel, pathOut = check_args(args)
    return opts, fileModel, pathOut

def check_args(args):
    if len(args) is not 2:
        usage("Improper number of arguments.")
    fileModel = args[0]
    pathOut = args[1]

    if not os.path.isfile(fileModel):
        usage("{0} is not a valid file".format(fileModel))

    pathsplit, dirsplit = os.path.split(pathOut)
    if not pathsplit:
        pathsplit = '.'
    if not os.path.isdir(pathsplit):
        usage("The path {0} does not exist".format(pathOut))

    if not os.path.isdir(pathOut):
        os.makedirs(pathOut)
    return fileModel, pathOut

# Print error messages and exit
def usage(errstr):
    print ""
    print "ERROR: {0}".format(errstr)
    print ""
    p.print_help()
    print ""
    exit(1)

def print_header():
    strDate = str(datetime.datetime.now())
    date, time = strDate.split()
    time = time.split('.')[0]
    print "                            _      ____  __  __ ______ " 
    print "                           | |    / __ \|  \/  |  ____|"
    print " _ __ ___   ___  _ __ _ __ | |__ | |  | | \  / | |__   "
    print "| '_ ` _ \ / _ \| '__| '_ \| '_ \| |  | | |\/| |  __|  "
    print "| | | | | | (_) | |  | |_) | | | | |__| | |  | | |____ "
    print "|_| |_| |_|\___/|_|  | .__/|_| |_|\____/|_|  |_|______|"
    print "                     | |                               "   
    print "                     |_|                               "
    print ""
    print "   Quantifying morphology from microscopic big data    "
    print "    https://github.com/slash-segmentation/morphome     "
    print ""
    print "                 {0} {1}".format(date, time)
    print ""

def nucleus(model, basename):
    scale = model.getScale()
    trans = model.getTrans()

    fnamewrl = os.path.join(pathOut, basename + '.wrl')
    fnamehx = os.path.join(pathOut, basename + '.hx')

    # Run nucleus pre-processing
    morphome.preprocess.nucleus(model, fnamewrl)

    # Copy nucleus workflow template .hx file to the output path. Replace
    # all instances of <FILENAME> with the appropriate file name.
    shutil.copyfile(os.path.join(pathHx, 'amira', 'proc_nucleus.hx'),
        fnamehx)


    regex_filename = re.compile(r"<FILENAME>")
    regex_pathout = re.compile(r"<PATH_OUT>")
    regex_scale = re.compile(r"<SCALE>")
    strScale = '{0} {1} {2}'.format(scale[0], scale[1], scale[2])
    fid2 = open(fnamehx + '.new', 'w')
    with open(fnamehx, 'rw') as fid:
        for line in fid:
            line = regex_filename.sub(fnamewrl, line)
            line = regex_pathout.sub(pathOut, line)
            line = regex_scale.sub(strScale, line)
            fid2.write(line)
    fid.close()
    fid2.close()
    shutil.move(fnamehx + '.new', fnamehx)

    # Run the script in Amira
    cmd = '{0} -no_gui {1}'.format(binAmira, fnamehx)
    subprocess.call(cmd.split())

    # Find nuclear envelope folds
    fshapeindex = os.path.join(pathOut, basename + '_shapeindex.am')
    fgrad = os.path.join(pathOut, basename + '_gradient.am')
    fsurf = os.path.join(pathOut, basename + '_surface.surf')
    fconvhull = os.path.join(pathOut, basename + '_convhull.surf')
    shapeIndex = morphome.readfile.scalar_field(fshapeindex)
    gradient = morphome.readfile.vector_field(fgrad)
    vertices, indices = morphome.readfile.surface(fsurf)
    vch, _ = morphome.readfile.surface(fconvhull)

    coords, sa = morphome.nucleus.quantify_nuclear_folds(model, vertices,
        indices, shapeIndex, gradient)

    # Do convex hull stuff
    nVertsCh = vch.shape[0]
    for i in range(nVertsCh):
        vch[i,:] = morphome.utils.transform_coords(model, vch[i,0], vch[i,1],
            vch[i,2])
        vch[i,:] = [int(x) for x in vch[i,:]] 
    vch = vch[vch[:,2].argsort()]

    for i in range(nVertsCh):
        print vch[i,:]


    pyimod.ImodWrite(model, 'blah.mod')

if __name__ == '__main__':
    global pathOut, pathHx, binAmira
    opts, fileModel, pathOut = parse_args()
    binAmira = '/usr/local/apps/Amira/6.0.0/bin/Amira'
    pathHx = os.path.split(morphome.__file__)[0]

    print_header()

    orgDict = {'nucleus': 'nucleus',
               'nuclei': 'nucleus',
               'nuc': 'nucleus'}

    # Load model file and print info
    print "Loading model file: {0}".format(fileModel)
    origModel = pyimod.ImodModel(fileModel)
    nObj = origModel.nObjects
    print "# Objects: {0}".format(nObj)
    print ""

    # Loop over all objects
    for i in range(nObj):
        namei = origModel.Objects[i].name
        print "Object {0}".format(i+1)
        print "============"
        print "Name: {0}".format(namei)
        print "# Contours: {0}".format(origModel.Objects[i].nContours)

        # If the current object name is empty, print a warning and skip it
        if not namei:
            print "WARNING: No object name provided. Skipping the object."
            continue

        # Check if the current object name is supported. If it is not, print a
        # warning and skip it. If it is, run the desired workflow.
        if orgDict.has_key(namei.lower()):
            workflow = orgDict[namei.lower()]
        else:
            print "WARNING: No existing workflow for object name {0}".format(namei)
            continue
        print "Running morphome workflow: {0}()".format(workflow)

        # Extract the object to a new model file
        modi = pyimod.ImodCmd(origModel, 'imodextract {0}'.format(i+1))

        # Construct basename
        basename = 'obj_' + str(i+1).zfill(4) + '_' + workflow

        # Run the corresponding function for the desired organelle workflow.
        try:
            func = globals()[workflow]
        except KeyError, e:
            print "WARNING: No function for {0} workflow. Skipping the object.".format(workflow)
        else:
            func(modi, basename)
