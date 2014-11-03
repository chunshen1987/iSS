#!/usr/bin/env python

import sys;
if sys.version[:3]=='1.5':
  from lib152 import *;
else:
  from lib import *;

from sys import argv;
from listR import stringToNumbers, biDifference, FLL, listFormCubicLatticeD, tracedSort, applyOrderList;

class control_variables: pass
control_variables.to_skip = {}; # This list is used to determine how many lines will be skipped at the beginning of the each data file. For example, if it is {1:1,2:3} meaning that 1 line in file1 will be skipped and 3 lines in file2 will be skipped.

def readFormatFile(filename):
  """ Read a format file and return a dictionary. The each line of the format file should be of the form:
    "variable name" @ file_id -> column_id, file_id -> column_id, ...
    It means that the "column_id"-th column(s) of file(s) with "file_id" corresponds to (components of) data with name "variable name". If no "->" is found, then all "file_id" are assumed to be 1. Those lines without "@" will be skipped.
    The return dictionary is of the format:
    {"variable name": [[file_id, column_id],...]}
  """

  # global control_variables

  D = {};
  for aLine in open(filename):
    if "@" not in aLine: # skip those lines without "@":
      print("vtk_converter::readFormatFile warning: no @ symbol found in the line:\n"+aLine+"\n and it will be skipped.\n");
      continue;
    [name, column_mapping] = [x.strip() for x in aLine.split("@")];
    if name=="skip":
      # this is a control variable
      control_variables.to_skip = dict([[int(y) for y in stringToNumbers(x,"->")] for x in column_mapping.split(",")]);
    else:
      D[name] = [[int(y) for y in stringToNumbers(x,"->")] for x in column_mapping.split(",")];
  return D;


def fmtToMapping(fmt):
  """ Return a dictionary of the format:
  {"t":f_time, "coord":f_coord, "scalar":{"scalar_name1":f_scalar1,...}, "vector":{"vector_name1":f_vector1,...}}
  where f_xxx are functions that when applied to the data read from files (see convertToVtk) give the corresponding values. "f_time", missing coordinates, and missing vector components will all return 0.
  Note that the code can be made more efficient by replacing all the occurence of fmt[xx] by a pre-evaluated variable.
  """
  D = {}; # the huge dictionary to be returned

  # function for time first:
  value_string = "0" if "t" not in fmt.keys() else "var["+str(fmt["t"][0][0]-1)+"]["+str(fmt["t"][0][1]-1)+"]"; # -1 because file_id specified starts from 1 instead of 0 as in python
  exec("f = lambda var: "+value_string);
  D["t"] = f;

  # coordinates:
  value_strings = ["0" if var not in fmt.keys() else "var["+str(fmt[var][0][0]-1)+"]["+str(fmt[var][0][1]-1)+"]" for var in ["x", "y", "z"]];
  exec("f = lambda var: ["+",".join(value_strings)+"]");
  D["coord"] = f;

  # get keys for scalar and vectors:
  all_keys = biDifference(fmt.keys(), ["x","y","z"]);
  scalar_keys = [];
  for aKey in all_keys:
    if len(fmt[aKey])==1: scalar_keys.append(aKey);
  vector_keys = biDifference(all_keys, scalar_keys);

  # scalars:
  D["scalar"] = {}; # a good start
  for aScalar in scalar_keys: # for one scalar, repeat what we did to "t":
    value_string = "var["+str(fmt[aScalar][0][0]-1)+"]["+str(fmt[aScalar][0][1]-1)+"]";
    exec("f = lambda var: "+value_string);
    D["scalar"][aScalar] = f;

  # vectors:
  D["vector"] = {}; # a good start
  for aVector in vector_keys: # for one vector, do the following:
    value_strings = ["0"]*3;
    for ii in range(len(fmt[aVector])):
      value_strings[ii] = "var["+str(fmt[aVector][ii][0]-1)+"]["+str(fmt[aVector][ii][1]-1)+"]";
    exec("f = lambda var: ["+",".join(value_strings)+"]");
    D["vector"][aVector] = f;

  return D;

def initialDataLists(D):
  """ Initialize grid, scalars, and vectors lists used in convertToVtk function. """
  grid = []; # simple
  scalars = [list(x) for x in [()]*len(D["scalar"].keys())]; # more complicated
  vectors = [list(x) for x in [()]*len(D["vector"].keys())]; # more complicated
  return [grid, scalars, vectors];

def fillData(grid, scalars, vectors, D, aLine):
  """ This function fill the data read from aLine to the lists grid, scalar, and vector, according to the regulated rule D. See convertToVtk for details. """
  grid.append(D["coord"](aLine));
  scalars_keys = D["scalar"].keys();
  for ii in range(len(scalars_keys)):
    scalars[ii].append(D["scalar"][scalars_keys[ii]](aLine));
  vectors_keys = D["vector"].keys();
  for ii in range(len(vectors_keys)):
    vectors[ii].append(D["vector"][vectors_keys[ii]](aLine));


def dumpVtk(filename, grid, scalars, vectors, D, structured="auto"):
  """ Write data to file. If structured="yes" then structured gird will be used; if structure="no" then unstructured grid will be used; if structured="auto" then structured grid will be used if "grid" form a cubic lattice (see listR.listFormCubicLattice for info on cubic lattice).
  """
  # setup grid
  structured = structured.lower();
  if structured == "yes": use_structured_grid = True;
  if structured == "no": use_structured_grid = False;
  if structured == "auto":
    form_lattice_D = listFormCubicLatticeD(grid);
    use_structured_grid = form_lattice_D["answer"];

  if use_structured_grid == False: # unstructured grid; simple
    grid_cmd = "UnstructuredGrid(grid)";
  else: # structured grid; more labor
    [grid, permutation] = tracedSort(grid); # sort it so x changes 1st, y changes 2nd, etc
    scalars = [applyOrderList(permutation, x) for x in scalars];
    vectors = [applyOrderList(permutation, x) for x in vectors];
    grid_cmd = "StructuredGrid("+str(form_lattice_D["dim"])+", "+str(grid)+")";

  # fill in objects
  object_strings = [];

  scalars_keys = D["scalar"].keys();
  for ii in range(len(scalars_keys)):
    object_strings.append("Scalars(scalars["+str(ii)+"], name='"+scalars_keys[ii]+"')");

  vectors_keys = D["vector"].keys();
  for ii in range(len(vectors_keys)):
    object_strings.append("Vectors(vectors["+str(ii)+"], name='"+vectors_keys[ii]+"')");

  # call VtkData
  exec("VtkData("+grid_cmd+", PointData("+','.join(object_strings)+")).tofile(filename)");


def convertToVtk(output_file, format_file, file_list, structured="auto"):
  """ Convert data given in the list "file_list" to "output_file". The name of "output_file" should contain a "%d" string for time slicing (t=0 if not included). Format of the files are defined in "format_file" (see readFormatFile function). The variable "x", "y", "z" are used to defined (unstructured) grids. The variable "t" is used to slice the data file: everytime its value changes, a new file "output_file % index" will be used for the output.
  """

  file_list = FLL([file_list]); # to make it more user friendly: for single file there need no "[]" around it
  fmt = readFormatFile(format_file); # raw format
  D = fmtToMapping(fmt); # regulated format
  t_id = 0;# controls the time-slice of the output

  # open files
  fids = [open(x) for x in file_list];
  # skip lines
  for skip_fid in control_variables.to_skip.keys():
    [fids[skip_fid-1].readline() for i in range(control_variables.to_skip[skip_fid])] # fid-1 b/c indices start from 0

  firstline_flag = 1;
  [grid, scalars, vectors] = initialDataLists(D); # grid is the collection of grid points, it is of the form [[x,y,z],...], the missed variables will be set to 0. scalars and vectors are nested lists, they has the format scalars=[[[#,#,...],[#,#,...],...],vectors=[[[#,#,#],[#,#,#],...],...]
  aStringLine = [x.readline() for x in fids];
  aLine = [stringToNumbers(y,",") for y in aStringLine]; # read one line from all files, then convert to nested lists. For example, if the 1st file has 2 columns and the 2nd has 1, then aLine will be like [[#,#],[#]]
  while aStringLine[0]!="": # total number of read rows is determined by the number of rows of the first file
    if aLine[0]!=[]: # not a comment-only line
      new_time=D["t"](aLine);
      if firstline_flag:
        old_time=new_time; # skip the 1st line
        firstline_flag=0; # no longer the 1st line
      if new_time!=old_time: # time step changed!
        dumpVtk(output_file % t_id, grid, scalars, vectors, D, structured); # output current slice
        old_time = new_time; # update current time
        [grid, scalars, vectors] = initialDataLists(D); # clear buffer
        fillData(grid, scalars, vectors, D, aLine); # add newly read line
        t_id=t_id+1; # advance to the next time step
      else: # time step kept; keep filling data
        fillData(grid, scalars, vectors, D, aLine);
    aStringLine = [x.readline() for x in fids];
    aLine = [stringToNumbers(y,",") for y in aStringLine]; # read one line from all files, then convert to nested lists. For example, if the 1st file has 2 columns and the 2nd has 1, then aLine will be like [[#,#],[#]]
  dumpVtk(output_file % t_id, grid, scalars, vectors, D, structured); # the final output
  print("Finished. Thanks for using vtk_converter. Zhi Qiu 2011");


if __name__ == '__main__':
  if len(argv)<4:
    print("Usage: vtk_converter.py output_file format_file file1 file2 ...");
  else:
    output_file = argv[1]+"_%d"; format_file = argv[2];
    convertToVtk(output_file, format_file, argv[3:], structured="auto");
