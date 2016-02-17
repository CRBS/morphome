#!/usr/bin/env python

import sys
import os
import pyimod
import morphome
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

def nucleus(model, basename):
    morphome.preprocess.nucleus(model, basename + '.wrl')

if __name__ == '__main__':
    opts, fileModel, pathOut = parse_args()

    print "============"
    print "= morphOME ="
    print "============"

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

