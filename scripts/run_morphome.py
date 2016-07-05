#!/usr/bin/env python

import sys
import os
import re
import glob
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
    p.add_option("--write_tiles_mitochondrion",
                 action = "store_true",
                 default = False,
                 dest = "writeTilesMitochondrion",
                 help = "Toggles on writing of tiled animations for individual "
                        "mitochondria.")
    p.add_option("--write_sql",
                 dest = "write_sql",
                 metavar = "FILE",
                 help = "Toggles on writing of SQLite output to the specified "
                        "file. If the file already exits, values will be added. "
                        "If not, a new database will be generated." )
    (opts, args) = p.parse_args()
    fileModel, pathOut = check_args(args)
    return opts, fileModel, pathOut

def check_args(args):
    if len(args) is not 2:
        usage("Improper number of arguments.")
    fileModel = args[0]
    pathOut = args[1]

    if not os.path.isfile(fileModel) and not os.path.isdir(fileModel):
        usage("{0} does not exist".format(fileModel))

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

def mitochondrion(model, basename, filename):
    scale = model.getScale()
    trans = model.getTrans()
 
    pathouti = os.path.join(pathOut, filename, basename)

    fnamewrl = os.path.join(pathouti, basename + '.wrl')
    fnamehx = os.path.join(pathouti, basename + '.hx')

    # Run mitochondrion pre-processing
    morphome.preprocess.mitochondrion(model, fnamewrl)

    shutil.copyfile(os.path.join(pathHx, 'amira', 'proc_mitochondrion.hx'),
        fnamehx)
    regex_filename = re.compile(r"<FILENAME>")
    regex_pathout = re.compile(r"<PATH_OUT>")
    regex_scale = re.compile(r"<SCALE>")
    regex_display = re.compile(r"<DISPLAY>")
    strScale = '{0} {1} {2}'.format(scale[0], scale[1], scale[2])
    fid2 = open(fnamehx + '.new', 'w')
    with open(fnamehx, 'rw') as fid:
        for line in fid:
            line = regex_filename.sub(fnamewrl, line)
            line = regex_pathout.sub(pathouti, line)
            line = regex_scale.sub(strScale, line)
            line = regex_display.sub(str(int(opts.writeTilesMitochondrion)), line)
            fid2.write(line)
    fid.close()
    fid2.close()
    shutil.move(fnamehx + '.new', fnamehx)

    # Run the script in Amira
    if opts.writeTilesMitochondrion:
        cmd = '{0} {1}'.format(binAmira, fnamehx)
    else:
        cmd = '{0} -no_gui {1}'.format(binAmira, fnamehx)
    subprocess.call(cmd.split())

def nucleus(model, basename, filename):
    scale = model.getScale()
    trans = model.getTrans()

    pathouti = os.path.join(pathOut, filename, basename)

    fnamewrl = os.path.join(pathouti, basename + '.wrl')
    fnamehx = os.path.join(pathouti, basename + '.hx')

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
            line = regex_pathout.sub(pathouti, line)
            line = regex_scale.sub(strScale, line)
            fid2.write(line)
    fid.close()
    fid2.close()
    shutil.move(fnamehx + '.new', fnamehx)

    # Run the script in Amira
    cmd = '{0} -no_gui {1}'.format(binAmira, fnamehx)
    subprocess.call(cmd.split())

    # Find nuclear envelope folds
    fshapeindex = os.path.join(pathouti, basename + '_shapeindex.am')
    fgrad = os.path.join(pathouti, basename + '_gradient.am')
    fsurf = os.path.join(pathouti, basename + '_surface.surf')
    fconvhull = os.path.join(pathouti, basename + '_convhull.surf')
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

def process_object(model, i, filename):
    """
    For an input model and object number, extracts the object to a new pyimod
    ImodModel class instance, determines the object name, and executes the
    appropriate workflow. If the object name is not part of the organelle
    dictionary, print a warning message and exit.   
    """
    namei = model.Objects[i].name
    print "Object {0}".format(i+1)
    print "============"
    print "Name: {0}".format(namei)
    print "# Contours: {0}".format(model.Objects[i].nContours)

    # If the current object name is empty, print a warning and skip it
    if not namei:
        print "WARNING: No object name provided. Skipping the object."
        return ''

    # Check if the current object name is supported. If it is not, print a
    # warning and skip it. If it is, run the desired workflow.
    if orgDict.has_key(namei.lower()):
        workflow = orgDict[namei.lower()]
    else:
        print "WARNING: No existing workflow for object name {0}".format(
            namei)
        return ''
    print "Running morphome workflow: {0}()".format(workflow)

    # Extract the object to a new, temporary model file
    modi = pyimod.ImodCmd(model, 'imodextract {0}'.format(i+1))

    # Construct basename to consist of the object name followed by the 
    # workflow/organelle name.
    basename = 'obj_' + str(i+1).zfill(4) + '_' + workflow

    # Make output directory for the object
    os.makedirs(os.path.join(pathOut, filename, basename))

    # Run the corresponding function for the desired organelle workflow.
    try:
        func = globals()[workflow]
    except KeyError, e:
        print "WARNING: No function for {0} workflow. Skipping the object.".format(workflow)
        return ''
    else:
        func(modi, basename, filename)
        return workflow

def write_csv_mitochondrion(filesIn):
    # Construct generic CSV header string
    headerstr = 'surface_area,' \
        'volume,' \
        'sav_ratio,' \
        'sphericity,' \
        'anisotropy,' \
        'centroid_x,' \
        'centroid_y,' \
        'centroid_z,' \
        'crofton_perimeter,' \
        'elongation,' \
        'euler_number_3d,' \
        'equivalent_diameter,' \
        'feret_shape_3d,' \
        'flatness,' \
        'shape_va3d,' \
        'orientation_phi,' \
        'orientation_theta,' \
        'n_branches,' \
        'total_length,' \
        'n_nodes,' \
        'n_intermediate_nodes,' \
        'n_terminal_nodes,' \
        'n_branching_nodes,' \
        'n_isolated_nodes' \

    # Construct a final header that is dependent on the maxmium number of
    # segments across all mitochondria computed. First, get a list of all
    # _length.csv files that contain branch data. Then, read each file and
    # extract the number of segments. Once all files have been read, take the
    # maximum number of segments across all mitochondria, and create a 
    # correspondingly sized string to append to the generic header string.    
    headerfiles = sorted(glob.glob(os.path.join(pathOut,
        '*/obj*mitochondrion/*_length.csv'))) 
    nseglist = []
    for fname in headerfiles:
        nSegments, _, _, _ = morphome.readfile.skel_length(fname)
        nseglist.append(nSegments)
    nsegmax = max(nseglist)
    for i in range(nsegmax):
        headeri = ',branch_{0}_length,' \
            'branch_{0}_mean_radius,' \
            'branch_{0}_orientation_theta,' \
            'branch_{0}_orientation_phi'.format(i+1)
        headerstr += headeri
    headerlen = headerstr.count(',') + 1

    for modelfile in filesIn:
        filecsv = os.path.join(pathOut, modelfile, 'mitochondrion.csv')
        fid = open(filecsv, 'a+')
        fid.write(headerstr + '\n') 
        objects = sorted(glob.glob(os.path.join(pathOut, modelfile, 'obj*mitochondrion')))
        for objectpath in objects:
            flength = glob.glob(os.path.join(objectpath, '*_length.csv'))[0]
            fsav = glob.glob(os.path.join(objectpath, '*_sav.csv'))[0]
            flabel = glob.glob(os.path.join(objectpath, '*_label.csv'))[0]
            nSegments, lengthTot, segments, nodes = morphome.readfile.skel_length(flength)
            sa, volume = morphome.readfile.sav(fsav)
            savratio, sphericity = morphome.utils.compute_sav_metrics(sa,
                volume)
            label_metrics = morphome.readfile.label_csv(flabel)

            # Convert metrics to microns as necessary
            for i in 1,2,3,4,7:
                label_metrics[i] = label_metrics[i] / 10000

            # Create final CSV file containing all computed metrics
            metrics = [sa, volume, savratio, sphericity]
            metrics.extend(label_metrics)
            metrics.extend([nSegments, lengthTot])
            metrics.extend(nodes)
            for i in range(nSegments):
                metrics.extend(segments[i,...])
            metrics.extend([None] * (headerlen - len(metrics)))
            fid.write(','.join([str(x) for x in metrics]) + '\n')
            os.remove(flength)
            os.remove(fsav)
            os.remove(flabel)
        fid.close()
    return filecsv

if __name__ == '__main__':
    global opts, pathOut, pathHx, binAmira, orgDict
    opts, fileModel, pathOut = parse_args()

    # Set path for Amira 6.0 on NCMIR machines
    binAmira = '/usr/local/apps/Amira/6.0.0/bin/Amira'
    pathHx = os.path.split(morphome.__file__)[0]

    # Print morphome header
    print_header()

    # Set dictionary to support singular/plural versions of organelle names
    # as well as a variety of common spellings.
    orgDict = {'nucleus': 'nucleus',
               'nuclei': 'nucleus',
               'mitochondrion': 'mitochondrion',
               'mitochondria': 'mitochondrion',
               'mito': 'mitochondrion'}

    # If the first supplied argument is a directory, parse the directory for
    # model files with the .mod extension. If it is a file, assume it is a
    # valid IMOD model file and process just this one file.
    if os.path.isdir(fileModel):
        filesIn = sorted(glob.glob(os.path.join(fileModel, '*.mod')))
    else:
        filesIn = [fileModel]
    
    # Loop over each model file. Get the number of objects in the file, and then
    # loop over each object, running the appropriate workflow.
    workflows_run = []
    for fname in filesIn:
        os.makedirs(os.path.join(pathOut, os.path.basename(fname)))
        print "Loading model file: {0}".format(fname)
        origModel = pyimod.ImodModel(fname)
        nObj = origModel.nObjects
        print "# Objects: {0}".format(nObj)
        print ""

        # Loop over each object and execute the appropriate workflow
        for iObj in range(nObj): 
            workflow_name = process_object(origModel, iObj, os.path.basename(
                fname))
            if workflow_name and (workflow_name not in workflows_run):
                workflows_run.append(workflow_name)

    for wname in workflows_run:
        try:
            func = globals()['write_csv_' + wname]
        except KeyError, e:
            print "WARNING: No function for {0} workflow. Skipping.".format(wname)
            break
        else:
            file_csv = func(filesIn)
            if opts.write_sql:
                print "Writing to SQL file {0}".format(opts.write_sql)
                conn, cur = morphome.sql.connect(opts.write_sql)
                cur = morphome.sql.insert_dataset(cur, (81739, 'Mouse', 'Suprachiasmatic Nucleus', 4, 1))
                cur = morphome.sql.insert_cell(cur, (2, 'Neuron'))
                #cur = morphome.sql.insert_organelle(cur, wname)
                cur = morphome.sql.read_csv(cur, wname, file_csv, 2)
                conn.commit()
                conn.close()
