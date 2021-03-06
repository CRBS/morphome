# Amira Project 600
# Amira
# Generated by morphOME (https://github.com/slash-segmentation/morphome)

proc makeCameraOrbitEvent {type time1 time2 angle1 angle2} {
    # Generates a single classic Demo Maker event string that performs rotation
    # of an object by a specified number of degrees across a given time 
    # interval.
    #   
    # Inputs
    # ------
    # type     String specifying the object type (i.e. "Camera-Orbit").
    # time1    Time to begin the rotation.
    # time2    Time to end the rotation (in seconds).
    # angle1   Angle to begin the rotation at (in degrees).
    # angle2   Angle to end the rotation at.
    #   
    # Returns
    # -------
    # str      Classic Demo Maker event string specifying the desired rotation.
    set str "dummy {numeric {$type/Time} $time1 $time2 $angle1 $angle2 0 360 "
    append str "{\"$type\" time setValue %0%; \"$type\" fire}} "
    return $str
}

proc calcShapeIndex {val1 val2} {
    # Smooths an input surface and computes its corresponding Shape Index
    # scalar field. Returns a value of 1 if the resultant scalar field does
    # not contain any NaNs. Otherwise, returns a zero.
    #
    # Inputs
    # ------
    # val1: Number of smoothing iterations
    # val2: Lambda value of smoothing
    #
    # Returns
    # -------
    # validSurf: = 0 if Shape Index scalar field has NaNs, else = 1.

    # Smooth the remeshed surface with the input parameters.
    set module "Smooth Surface"
    create HxSurfaceSmooth $module
    $module data connect "GeometrySurface.remeshed"
    $module parameters setValue 0 $val1
    $module parameters setValue 1 $val2
    $module action snap
    $module fire

    # Calculate the Shape Index on each vertex of the smoothed surface.
    set module "Curvature"
    create HxGetCurvature $module
    $module data connect "GeometrySurface.smooth"
    $module method setValue 1
    $module output setValue 10  
    $module doIt snap
    $module fire

    # Check that the Shape Index scalar field does not contain any NaNs. If it
    # does, this is typically indicative of an overly smoothed surface, and the
    # smoothing will need to be re-run with more conservative parameters.
    set validSurf 1
    set nNodes ["ShapeIndex" nValues]
    for {set i 0} {$i < $nNodes} {incr i} {
        set valSI ["ShapeIndex" getValue $i] 
        if {[isnan $valSI]} {
            set validSurf 0
            remove "ShapeIndex" "Curvature"
            remove "GeometrySurface.smooth" "Smooth Surface"
            break
        }
    }   
    return $validSurf
}

proc isnan {x} {
    # Checks if the input is NaN.
    #
    # Input
    # -----
    # x: Variable to be tested
    #
    # Returns
    # -------
    # = 1 if x is NaN, else = 0
    
    if {![string is double $x] | $x != $x} {
        return 1
    } else {
        return 0
    }   
}

###
#####
# HARD-CODED INPUT VALUES
#####
###

set fname "<FILENAME>"
set pathOut "<PATH_OUT>"
set displayOn "<DISPLAY>"
set scale [list <SCALE>]
set scalex [lindex $scale 0]
set scaley [lindex $scale 1]
set scalez [lindex $scale 2]

###
#####
# DATA IMPORT
#####
###

# Extract the VRML file's basename
set base [file tail $fname]
set base [string trimright $base ".wrl"]

# Load the VRML file into Amira
[load $fname] setLabel $base

# Convert VRML to surface (named "GeometrySurface")
set module "Open Inventor Scene To Surface"
create HxGeometryToSurface $module
$module data connect $base
$module action snap
$module fire

# Remesh the surface. The surface is remeshed using the best isotropic vertex
# placement and 3x the number of input triangles. This is necessary to get a 
# finer sampling of the Shape Index across the surface.
set module "Remesh Surface"
set nTriIn ["GeometrySurface" getNumTriangles]
set nTriOut [expr $nTriIn * <SURFACE_TRIANGLE_MULTIPLIER>]
create HxRemeshSurface $module
$module select
$module data connect "GeometrySurface"
$module fire
$module objective setIndex 0 1 
$module interpolateOrigSurface setValue 0
$module desiredSize setValue 1 $nTriOut
$module remeshOptions1 setValue 0 0 
$module remeshOptions1 setValue 1 1 
$module fire
$module remesh snap
$module fire

###
#####
# SHAPE INDEX & NUCLEAR FOLDING
#####
###

# Smooth the remeshed surface and compute the Shape Index. If the smoothed
# surface is invalid, retry smoothing with less liberal smoothing parameters.
# As a first attempt, the number of iterations will be decreased. If this is
# unsuccessful, the lambda value will be decreased. If no valid smoothed
# surface can be found, the program will print a warning and exit.
set validSurf 0
set smoothVal1 <SURFACE_SMOOTH_ITERATIONS>
set smoothVal2 <SURFACE_SMOOTH_LAMBDA>
while {!$validSurf} {
    echo "Trying smoothing with $smoothVal1, $smoothVal2"
    set validSurf [calcShapeIndex $smoothVal1 $smoothVal2]
    if {!$validSurf} {
        if {$smoothVal1 > 2} {
            set smoothVal1 [expr $smoothVal1 - 1]
        } elseif {$smoothVal1 <= 2 & $smoothVal2 > 0.2} {
            set smoothVal2 [expr $smoothVal2 - 0.1]
        } else {
            echo "WARNING: Smoothing failed. Exiting"
            exit
        }
    }   
}

# Compute gradient of the Shape Index scalar field.
set module "Surface Gradient"
create HxComputeSurfaceGradient $module
$module data connect "ShapeIndex"
$module fire

# Convert the gradient vectors from Angstroms to microns
set module "Arithmetic"
create HxArithmetic $module
$module inputA connect "ShapeIndex_Gradient"
$module expr0 setValue "Ax * 10000"
$module expr1 setValue "Ay * 10000"
$module expr2 setValue "Az * 10000"
$module create

###
#####
# WILLMORE ENERGY
#####
###

# Compute mean curvature on triangle faces
set module "Curvature 2"
create HxGetCurvature $module
$module data connect "GeometrySurface.smooth"
$module method setValue 0
$module output setValue 2
$module doIt snap
$module fire

# Compute Gaussian curvature on triangle faces
set module "Curvature 3"
create HxGetCurvature $module
$module data connect "GeometrySurface.smooth"
$module method setValue 0
$module output setValue 4
$module doIt snap
$module fire

###
#####
# CONVEX HULL DIFFERENCE
##### 
###

# Convex hull
set module "Convex Hull"
create HxConvexHull $module
"Convex Hull" data connect "GeometrySurface.smooth"
"Convex Hull" action snap
"Convex Hull" fire

# Get surface area and volume of nucleus
set module "Surface Area Volume Nucleus"
create HxSurfaceArea $module
$module data connect "GeometrySurface.smooth"
$module doIt snap
$module fire
set saNucleus ["GeometrySurface.statistics" getValue 2 0]
set vNucleus ["GeometrySurface.statistics" getValue 3 0]

# Get surface area and volume of convex hull
set module "Surface Area Volume Convex Hull"
create HxSurfaceArea $module
$module data connect "GeometrySurface-convexHull"
$module doIt snap
$module fire
set saConvHull ["GeometrySurface-convexHull.statistics" getValue 2 0]
set vConvHull ["GeometrySurface-convexHull.statistics" getValue 3 0]

###
#####
# 3D BINARY IMAGE METRICS
#####
###

# Convert the surface to a binary volume stack. First, the dimensions of the
# output stack need to be set. If the default dimensions aren't changed, Amira
# will try to output a stack requiring petabytes of memory.
set module "Scan Surface To Volume"
create HxScanConvertSurface $module
$module data connect "GeometrySurface.smooth"
$module field disconnect
$module fire
set xmin [$module bbox getValue 0]
set xmax [$module bbox getValue 1]
set ymin [$module bbox getValue 2]
set ymax [$module bbox getValue 3]
set zmin [$module bbox getValue 4]
set zmax [$module bbox getValue 5]
set dimx [expr round((double($xmax) / $scalex - double($xmin) / $scalex))]
set dimy [expr round((double($ymax) / $scaley - double($ymin) / $scaley))]
set dimz [expr round((double($zmax) / $scalez - double($zmin) / $scalez))]
$module dimensions setValues $dimx $dimy $dimz
$module action snap
$module fire
"GeometrySurface.scanConverted" sharedColormap setValue "grayScale.am"
"GeometrySurface.scanConverted" sharedColormap setMinMax 0 1
"GeometrySurface.scanConverted" ImageData disconnect
"GeometrySurface.scanConverted" fire
"GeometrySurface.scanConverted" primary setIndex 0 0
"GeometrySurface.scanConverted" fire
"GeometrySurface.scanConverted" select
"GeometrySurface.scanConverted" VoxelSize setValue "$scalex x $scaley x $scalez"

# Label analysis
set module "Label Analysis"
create HxAnalyzeLabels $module
$module data connect "GeometrySurface.scanConverted"
$module measures setState "Nucleus" Anisotropy BaryCenterX BaryCenterY \
    BaryCenterZ CroftonPerimeter Elongation Euler3D EqDiameter FeretShape3d \
    Flatness Shape_VA3d OrientationPhi OrientationTheta
$module interpretation setValue 0
$module doIt hit 
$module fire

###
#####
## DATA EXPORT
#####
###

## Export files to disk for further analysis
set fname_gradient [file join $pathOut ${base}_gradient.am]
set fname_shapeindex [file join $pathOut ${base}_shapeindex.am]
set fname_surface [file join $pathOut ${base}_surface.surf]
set fname_labelcsv [file join $pathOut ${base}_label.csv]
set fname_savcsv_nuc [file join $pathOut ${base}_sav_nuc.csv]
set fname_savcsv_ch [file join $pathOut ${base}_sav_ch.csv]
set fname_meancurv [file join $pathOut ${base}_meancurv.am]
set fname_gausscurv [file join $pathOut ${base}_gausscurv.am]

"Result" exportData "Amira ASCII" $fname_gradient
"ShapeIndex" exportData "Amira ASCII" $fname_shapeindex
"MeanCurvature" exportData "Amira ASCII" $fname_meancurv
"GaussCurvature" exportData "Amira ASCII" $fname_gausscurv
"GeometrySurface.smooth" exportData "HxSurface ASCII" $fname_surface
"GeometrySurface.Label-Analysis" exportData "CSV" $fname_labelcsv
"GeometrySurface.statistics" exportData "CSV" $fname_savcsv_nuc
"GeometrySurface-convexHull.statistics" exportData "CSV" $fname_savcsv_ch

###
#####
## DISPLAY
#####
##

if {$displayOn} {
    # Set background to solid black
    viewer 0 setBackgroundMode 0
    viewer 0 setBackgroundColor 0 0 0

    # Set up display of the surface
    set module "Surface View"
    create HxDisplaySurface $module
    $module data connect "GeometrySurface.smooth"
    $module drawStyle setValue 4
    $module drawStyle setNormalBinding 1
    $module colorMode setValue 5
    $module colormap setDefaultColor <SURFACE_RED> <SURFACE_GREEN> <SURFACE_BLUE>
    $module baseTrans setValue 0
    $module fire

    # Create camera path, rotating once around the object
    set module "Camera-Orbit"
    create HxCircularCameraPath $module
    $module action setState menus 1 4 
    $module fire

    # Setup Demo Maker
    set movieLength 4
    set movieEventString [makeCameraOrbitEvent "Camera-Orbit" 0 4 0 360]
    echo $movieEventString
    set scriptAnim [load ${AMIRA_ROOT}/share/script-objects/DemoMakerClassic.scro]
    set scriptAnim [lindex $module 0]
    $scriptAnim setVar scroTypeDemoMaker 1
    $scriptAnim setVar "internalEventList" $movieEventString
    $scriptAnim setVar "lastStartTime" 0
    $scriptAnim setVar "lastEndTime" $movieLength
    $scriptAnim setVar "lastTimeStep" 0
    $scriptAnim setVar "loadNetwDemoMakers" {}
    $scriptAnim fire
    $scriptAnim time disconnect
    $scriptAnim time setMinMax 0 $movieLength
    $scriptAnim time setSubMinMax 0 $movieLength
    $scriptAnim time setValue 1
    $scriptAnim time setValue 0
    $scriptAnim time animationMode -once
    $scriptAnim fire
    $scriptAnim select

    # Create Movie Maker module
    set module "Movie Maker"
    create HxMovieMaker $module
    $module time connect $scriptAnim
    $module fire
    $module fileFormat setValue 2
    $module frameRate setValue 0
    $module frames setValue 18
    $module compressionQuality setValue 1
    $module size setValue 5
    $module resolution setValue 0 500
    $module resolution setValue 1 500
    $module filename setValue $pathOut/mov
    $module action setState index 0
    $module action touch 0
    $module fire
    echo $fname
    echo $pathOut
}

exit
