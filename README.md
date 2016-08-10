[numpy]: http://www.numpy.org/
[pyimod]: https://github.com/CRBS/PyIMOD
[amira]: https://www.fei.com/software/amira-3d-for-life-sciences/
[imod]:http://bio3d.colorado.edu/imod/
[vrml]:https://en.wikipedia.org/wiki/VRML

# morphome
**morphome** is a set of Python and Amira Tcl scripts designed to easily and automatically import [IMOD][imod] 3D model files into FEI Amira and perform organelle-specific morphological quantification routines. **morphome** uses the [PyIMOD][pyimod] set of Python classes to read an input IMOD model file and convert its objects into individual 3D files in the Virtual Reality Modeling Language ([VRML][vrml]) format. Once conversion is complete, each individual VRML file is read into Amira and individual quantification routines are launced depending on the object name specified in the IMOD model file.

Currently, **morphome** supports quantification workflows for the following organelles:
* Mitochondrion
* Nucleus
* Nucleolus

### Requirements
* Python 2.7
* [Numpy][numpy]
* [PyIMOD][pyimod]
* [IMOD][imod]
* [Amira 6.0][amira]

### Example Usage
To run using the default settings:  

    morphome/scripts/run_morphome.py input_file.mod /home/user/output_path

To run using custom settings stored to a JSON file:

    morphome/scripts/run_morphome.py --json custom_settings.json input_file.mod /home/user/output_path
