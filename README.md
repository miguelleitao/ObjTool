![image](https://travis-ci.org/miguelleitao/pipeline.svg?branch=master "Linux Build Status") 
# ObjTool
The Swiss Army Knife for Wavefront OBJ models

ObjTool is a command-line utility to work with Wavefront OBJ files.

Uses OBJ file format as defined in http://paulbourke.net/dataformats/obj/
Not all OBJ file format specs are implemented.

_________
### Features
* Join geometry from several OBJ files into one.
* Select (use only) geometry that is (or is not) part of a specific object.
* Select (use only) geometry that is (or is not) marked with a specific material.
* Select (use only) geometry according to its position.
* Apply transformations (translation/rotation/scale) to selected (or all) geometry.

_________
### Dependences
**ObjTool** is implemented in C from scratch. It does not have any special dependences.
To build ObjTool you will only require a C compiler (usually gcc) and the linmath library.

_________
### Download / Build / Install

    git clone --recurse-submodules git@github.com:miguelleitao/ObjTool.git
    cd ObjTool
    make
    sudo make install
    ObjTool

_________
### Usage
```
ObjTool [options] inputfile1.obj inputfile2.obj
    -i                  file info
    -xmin value         define xmin
    -xmax value         define xmax
    -ymin value         define ymin
    -ymax value         define ymax
    -zmin value         define zmin
    -zmax value         define zmax
    -o name             select object
    -g name 	        select obj group
    -m name             select material
    -tx value	        translate x
    -ty value	        translate y
    -tz value	        translate z
    -sx value	        scale x
    -sy value	        scale y
    -sz value	        scale z
    -s value	        global scale
    -rx value	        rotate around xx axis
    -ry value	        rotate around yy axis
    -rz value	        rotate around zz axis
    -N			invert normals
    -R                  use relative coords    
    -M	                Do not output mtllib directives
    -c                  solid cut
    -n                  negate face filter condition( -g,-m )
    -O outfile          output to outfile (default: stdout)
    -e                  Explode outfile into single objects
    -S shadow_file      shadow output to shadow_file (default: no shadow ouput)
```
