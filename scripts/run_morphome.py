#!/usr/bin/env python

import sys
import os
import re
import glob
import shutil
import subprocess
import datetime
import json
import numpy as np
from optparse import OptionParser

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import morphome
import pyimod

def parse_args():
    global p
    p = OptionParser(usage = "%prog [options] file.mod /path/for/output")

    p.add_option("--json",
                 dest = "json",
                 default = os.path.join(os.path.split(morphome.__file__)[0],
                     'json', 'defaults.json'),
                 metavar = "FILE",
                 help = "JSON file to read default values from.") 
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
    print "           https://github.com/crbs/morphome     "
    print ""
    print "                 {0} {1}".format(date, time)
    print ""

def mitochondrion(model, basename, filename):
    scale = model.getScale()
    trans = model.getTrans()

    # Get absolute paths for output VRML and .hx files  
    pathouti = os.path.join(pathOut, filename, basename)
    fnamewrl = os.path.abspath(os.path.join(pathouti, basename + '.wrl'))
    fnamehx = os.path.abspath(os.path.join(pathouti, basename + '.hx'))

    # Run mitochondrion pre-processing and convert to VRML
    morphome.preprocess.mitochondrion(model, fnamewrl)

    # Copy mitochondrion .hx template file from the morphome/amira directory
    # to the output path for the object's .hx file
    shutil.copyfile(os.path.join(pathHx, 'amira', 'proc_mitochondrion.hx'),
        fnamehx)

    # Compile regular expressions to match within the Amira script
    regex_filename = re.compile(r"<FILENAME>")
    regex_pathout = re.compile(r"<PATH_OUT>")
    regex_scale = re.compile(r"<SCALE>")
    regex_display = re.compile(r"<DISPLAY>")
    regex_surfr = re.compile(r"<SURFACE_RED>")
    regex_surfg = re.compile(r"<SURFACE_GREEN>")
    regex_surfb = re.compile(r"<SURFACE_BLUE>")
    regex_surft = re.compile(r"<SURFACE_TRIANGLE_MULTIPLIER>")
    regex_surfsi = re.compile(r"<SURFACE_SMOOTH_ITERATIONS>")
    regex_surfsl = re.compile(r"<SURFACE_SMOOTH_LAMBDA>")
    regex_skelw = re.compile(r"<SKELETON_WIDTH>")
    regex_skelc1 = re.compile(r"<SKELETON_SMOOTH_COEFFICIENT1>")
    regex_skelc2 = re.compile(r"<SKELETON_SMOOTH_COEFFICIENT2>")
    regex_skelsi = re.compile(r"<SKELETON_SMOOTH_ITERATIONS>")
    regex_noder = re.compile(r"<NODE_RED>")
    regex_nodeg = re.compile(r"<NODE_GREEN>")
    regex_nodeb = re.compile(r"<NODE_BLUE>")

    # Parse through Amira script, replacing regular expression matches with the
    # appropriate data supplied by the loaded JSON file.
    strScale = '{0} {1} {2}'.format(scale[0], scale[1], scale[2])
    data = json_data["mitochondrion"]    

    # Run through the copied .hx file line-by-line, and replace regular
    # expression matches when encountered.
    fid2 = open(fnamehx + '.new', 'w')
    with open(fnamehx, 'rw') as fid:
        for line in fid:
            line = regex_filename.sub(fnamewrl, line)
            line = regex_pathout.sub(os.path.abspath(pathouti), line)
            line = regex_scale.sub(strScale, line)
            line = regex_display.sub(str(data["write_animation"]), line)
            line = regex_surfr.sub(str(data["surface_red"]), line)
            line = regex_surfg.sub(str(data["surface_green"]), line)
            line = regex_surfb.sub(str(data["surface_blue"]), line)
            line = regex_surft.sub(str(data["surface_triangle_multiplier"]),
                line)
            line = regex_surfsi.sub(str(data["surface_smooth_iterations"]),
                line)
            line = regex_surfsl.sub(str(data["surface_smooth_lambda"]), line)
            line = regex_skelw.sub(str(data["skeleton_width"]), line)
            line = regex_skelc1.sub(str(data["skeleton_smooth_coefficient1"]),
                line)
            line = regex_skelc2.sub(str(data["skeleton_smooth_coefficient2"]),
                line)
            line = regex_skelsi.sub(str(data["skeleton_smooth_iterations"]),
                line)
            line = regex_noder.sub(str(data["node_red"]), line)
            line = regex_nodeg.sub(str(data["node_green"]), line)
            line = regex_nodeb.sub(str(data["node_blue"]), line) 
            fid2.write(line)
    fid.close()
    fid2.close()
    shutil.move(fnamehx + '.new', fnamehx)

    # Run the script in Amira
    if data["write_animation"]:
        cmd = '{0} {1}'.format(binAmira, fnamehx)
    else:
        cmd = '{0} -no_gui {1}'.format(binAmira, fnamehx)
    subprocess.call(cmd.split())

    # Remove the last 'exit' line of the .hx script
    remove_last_line(fnamehx)

def nucleus(model, basename, filename):
    scale = model.getScale()
    trans = model.getTrans()

    # Get absolute paths for output VRML, .hx, and IMOD model files
    pathouti = os.path.abspath(os.path.join(pathOut, filename, basename))
    fnamewrl = os.path.abspath(os.path.join(pathouti, basename + '.wrl'))
    fnamehx = os.path.abspath(os.path.join(pathouti, basename + '.hx'))
    fnamemod = os.path.abspath(os.path.join(pathouti, basename + '.mod'))

    # Run nucleus pre-processing and convert to VRML
    morphome.preprocess.nucleus(model, fnamewrl)

    # Copy nucleus workflow template .hx file to the output path. Replace
    # all instances of <FILENAME> with the appropriate file name.
    shutil.copyfile(os.path.join(pathHx, 'amira', 'proc_nucleus.hx'),
        fnamehx)

    # Compile regular expressions to match within the Amira script
    regex_filename = re.compile(r"<FILENAME>")
    regex_pathout = re.compile(r"<PATH_OUT>")
    regex_scale = re.compile(r"<SCALE>")
    regex_display = re.compile(r"<DISPLAY>")
    regex_surfr = re.compile(r"<SURFACE_RED>")
    regex_surfg = re.compile(r"<SURFACE_GREEN>")
    regex_surfb = re.compile(r"<SURFACE_BLUE>")
    regex_surft = re.compile(r"<SURFACE_TRIANGLE_MULTIPLIER>")
    regex_surfsi = re.compile(r"<SURFACE_SMOOTH_ITERATIONS>")
    regex_surfsl = re.compile(r"<SURFACE_SMOOTH_LAMBDA>")

    # Parse through Amira script, replacing regular expression matches with the
    # appropriate data supplied by the loaded JSON file.
    strScale = '{0} {1} {2}'.format(scale[0], scale[1], scale[2])
    data = json_data["nucleus"]

    # Run through the copied .hx file line-by-line, and replace regular
    # expression matches when encountered.
    fid2 = open(fnamehx + '.new', 'w')
    with open(fnamehx, 'rw') as fid:
        for line in fid:
            line = regex_filename.sub(fnamewrl, line)
            line = regex_pathout.sub(pathouti, line)
            line = regex_scale.sub(strScale, line)
            line = regex_display.sub(str(data["write_animation"]), line)
            line = regex_surfr.sub(str(data["surface_red"]), line)
            line = regex_surfg.sub(str(data["surface_green"]), line)
            line = regex_surfb.sub(str(data["surface_blue"]), line)
            line = regex_surft.sub(str(data["surface_triangle_multiplier"]),
                line)
            line = regex_surfsi.sub(str(data["surface_smooth_iterations"]),
                line)
            line = regex_surfsl.sub(str(data["surface_smooth_lambda"]), line)
            fid2.write(line)
    fid.close()
    fid2.close()
    shutil.move(fnamehx + '.new', fnamehx)

    # Run the script in Amira
    if data["write_animation"]:
        cmd = '{0} {1}'.format(binAmira, fnamehx)
    else:
        cmd = '{0} -no_gui {1}'.format(binAmira, fnamehx)
    subprocess.call(cmd.split())

    # Remove the last 'exit' line of the .hx script
    remove_last_line(fnamehx)

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
        return modi, workflow

def get_generic_header():
    """
    Returns a generic CSV header string that will be used for all organelle
    types.
    """ 
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
        'orientation_theta,'
    return headerstr

def get_generic_metrics(fsav, flabel):
    """
    Reads surface area/volume and label files output from Amira and returns
    metrics from each. The metrics read here are computed for every organelle
    workflow, and so are generic across all. 
    """
    # Read SA and volume data 
    sa, volume = morphome.readfile.sav(fsav)

    # Compute metrics from SA and volume data 
    savratio, sphericity = morphome.utils.compute_sav_metrics(sa,
        volume)

    # Read Amira label metrics
    label_metrics = morphome.readfile.label_csv(flabel)

    # Read Amira label metrics label_metrics = morphome.readfile.label_csv(flabel) 
    # Convert metrics to microns as necessary
    for i in 1,2,3,4,7:
        label_metrics[i] = label_metrics[i] / 10000

    return sa, volume, savratio, sphericity, label_metrics

def write_csv_mitochondrion(model, filesIn):
    """
    Compiles an output CSV file for mitochondria metrics across all mito objects
    and model files encountered. 
    """

    # Construct generic CSV header string consisting of the metrics that were
    # computed for every individual mitochondrion. 
    headerstr = get_generic_header()
    headerstr += 'n_branches,' \
        'total_length,' \
        'n_nodes,' \
        'n_intermediate_nodes,' \
        'n_terminal_nodes,' \
        'n_branching_nodes,' \
        'n_isolated_nodes' \

    # Construct a final header that is dependent on the maxmium number of
    # segments across all mitochondria computed. First, get a list of all
    # *_length.csv files that contain branch data. Then, read each file and
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

    # Loop over all files supplied. For each file, create a new CSV file
    # specific to mitochondria, and store it within the file's output path.
    # For each file, parse Amira outputs for each mito object, and append to
    # the CSV file.  
    for modelfile in filesIn:
        # Create a new CSV file for each file supplied, and open it to append to
        modelname = os.path.split(modelfile)[-1] 
        filecsv = os.path.join(pathOut, modelname, 'mitochondrion.csv')
        print "Writing mitochondrion metrics for {0} to {1}".format(modelfile,
            filecsv)
        fid = open(filecsv, 'a+')

        # Write global CSV header to the file 
        fid.write(headerstr + '\n') 

        # Find all mito objects for the given model file. Loop over each object.
        objects = sorted(glob.glob(os.path.join(pathOut, modelname,
            'obj*mitochondrion')))
        for objectpath in objects:
            # Get Amira output CSV file names
            flength = glob.glob(os.path.join(objectpath, '*_length.csv'))[0]
            fsav = glob.glob(os.path.join(objectpath, '*_sav.csv'))[0]
            flabel = glob.glob(os.path.join(objectpath, '*_label.csv'))[0]

            # Get generic metrics
            sa, volume, savratio, sphericity, label_metrics = \
                get_generic_metrics(fsav, flabel)

            # Read skeleton length data
            nSegments, lengthTot, segments, nodes = morphome.readfile.skel_length(flength)

            # Create a list containing all computed metrics
            metrics = [sa, volume, savratio, sphericity]
            metrics.extend(label_metrics)
            metrics.extend([nSegments, lengthTot])
            metrics.extend(nodes)
            for i in range(nSegments):
                metrics.extend(segments[i,...])
            metrics.extend([None] * (headerlen - len(metrics)))

            # Append the metrics list to the growing CSV file
            fid.write(','.join([str(x) for x in metrics]) + '\n')

            # Remove intermediate CSV files output by Amira
            os.remove(flength)
            os.remove(fsav)
            os.remove(flabel)
        fid.close()
    return filecsv

def write_csv_nucleus(model, filesIn):
    """ 
    Compiles an output CSV file for nucleus metrics across all mito objects and
    model files encountered. 
    """

    # Construct generic CSV header string consisting of the metrics that were
    # computed for every individual nuclei.
    headerstr = get_generic_header()
    headerstr += 'n_folds,' \
        'surface_area_folds,' \
        'fold_area_ratio,' \
        'volume_convex_hull,' \
        'convex_hull_difference,' \
        'willmore_energy'

    for modelfile in filesIn:
        # Create a new CSV file for each file supplied, and open it to append to
        modelname = os.path.split(modelfile)[-1] 
        filecsv = os.path.join(pathOut, modelname, 'nucleus.csv')
        print "Writing nucleus metrics for {0} to {1}".format(modelfile,
            filecsv)
        fid = open(filecsv, 'a+')

        # Write CSV header to the file 
        fid.write(headerstr + '\n') 

        # Find all nucleus objects for the given model file. Loop over each 
        # object.
        objects = sorted(glob.glob(os.path.join(pathOut, modelname,
            'obj*nucleus')))
        for objectpath in objects:
            basename = os.path.split(objectpath)[-1]

            # Get Amira output file names
            fsav = glob.glob(os.path.join(objectpath, '*_sav_nuc.csv'))[0]
            flabel = glob.glob(os.path.join(objectpath, '*_label.csv'))[0]
            fshapeindex = glob.glob(os.path.join(objectpath,
                '*_shapeindex.am'))[0]
            fgrad = glob.glob(os.path.join(objectpath, '*_gradient.am'))[0]
            fsurf = glob.glob(os.path.join(objectpath, '*_surface.surf'))[0]
            fconvhull = glob.glob(os.path.join(objectpath,
                '*_sav_ch.csv'))[0]
            fgauss = glob.glob(os.path.join(objectpath, '*_gausscurv.am'))[0]
            fmean = glob.glob(os.path.join(objectpath, '*_meancurv.am'))[0]

            # Get generic metrics
            sa, volume, savratio, sphericity, label_metrics = \
                get_generic_metrics(fsav, flabel)

            # Read relevant data from shape index and gradient scalar fields,
            # as well as mesh vertex and index data 
            shapeIndex = morphome.readfile.scalar_field(fshapeindex)
            gradient = morphome.readfile.vector_field(fgrad)
            vertices, indices = morphome.readfile.surface(fsurf)

            # Identify the nuclear envelope folds. Returns the number of folds
            # and surface area of folded regions of the nuclear surface. The 
            # ratio of folded to unfolded area is then computed.
            n_folds, sa_folds = morphome.nucleus.quantify_nuclear_folds(model,
                vertices, indices, shapeIndex, gradient)
            fold_sa_ratio = sa_folds / sa

            # Write IMOD model file with points added to denote fold locations.
            fnamemod = os.path.join(objectpath, basename + '.mod')
            pyimod.ImodWrite(model, fnamemod)

            # Read volume data from convex hull file. Compute convex hull 
            # difference.
            _, vol_ch = morphome.readfile.sav(fconvhull)
            ch_diff = vol_ch - volume 

            # Compute Willmore energy
            willmore_energy = morphome.nucleus.willmore_energy(fmean, fgauss,
                fsurf)       

            # Create a list containing all computed metrics
            metrics = [sa, volume, savratio, sphericity]
            metrics.extend(label_metrics)
            metrics.extend([n_folds, sa_folds, fold_sa_ratio, vol_ch, ch_diff,
                willmore_energy])

            # Append the metrics list to the growing CSV file
            fid.write(','.join([str(x) for x in metrics]) + '\n')

            # Remove intermediate CSV files output by Amira
            os.remove(fsav)
            os.remove(flabel)
            os.remove(fshapeindex)
            os.remove(fgrad)
            os.remove(fsurf)
            os.remove(fconvhull)
            os.remove(fgauss)
            os.remove(fmean)
        fid.close()
    return filecsv

def remove_last_line(fname):
    """
    Removes the last line of a supplied text file. In this case, the file is an
    Amira .hx script and the final line is an 'exit' command. This will be 
    removed so that a user can open the script directly after morphome is done
    running, and visualize the results in Amira.  
    """
    fid = open(fname, 'r+')
    fid.seek(0, os.SEEK_END)
    pos = fid.tell() - 1
    while pos > 0 and fid.read(1) != "\n":
        pos -= 1
        fid.seek(pos, os.SEEK_SET)
    if pos > 0:
        fid.seek(pos, os.SEEK_SET)
        fid.truncate()
    fid.close() 

if __name__ == '__main__':
    global opts, pathOut, pathHx, binAmira, orgDict, json_data
    opts, fileModel, pathOut = parse_args()

    # Set path for Amira 6.0 on NCMIR machines
    binAmira = '/usr/local/apps/Amira/6.0.0/bin/Amira'
    pathHx = os.path.split(morphome.__file__)[0]
    print pathHx

    # Print morphome header
    print_header()

    # Set dictionary to support singular/plural versions of organelle names
    # as well as a variety of common spellings.
    orgDict = {'nucleus': 'nucleus',
               'nuclei': 'nucleus',
               'mitochondrion': 'mitochondrion',
               'mitochondria': 'mitochondrion',
               'mito': 'mitochondrion'}

    # Read JSON settings file
    print "Reading JSON settings file: {}".format(opts.json)
    with open(opts.json) as file_json:
        json_data = json.load(file_json)   

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
            modi, workflow_name = process_object(origModel, iObj, 
                os.path.basename(fname))
            if workflow_name and (workflow_name not in workflows_run):
                workflows_run.append(workflow_name)

    for wname in workflows_run:
        try:
            func = globals()['write_csv_' + wname]
        except KeyError, e:
            print "WARNING: No function for {0} workflow. Skipping.".format(wname)
            break
        else:
            file_csv = func(modi, filesIn)
            if opts.write_sql:
                print "Writing to SQL file {0}".format(opts.write_sql)
                conn, cur = morphome.sql.connect(opts.write_sql)
                cur = morphome.sql.insert_dataset(cur, (81739, 'Mouse', 'Suprachiasmatic Nucleus', 4, 1))
                cur = morphome.sql.insert_cell(cur, (2, 'Neuron'))
                cur = morphome.sql.read_csv(cur, wname, file_csv, 2)
                conn.commit()
                conn.close()
