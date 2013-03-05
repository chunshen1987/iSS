Last modified: 08-28-2011, Zhi Qiu

This program converts a column-based data file (i.e. each column represents one type of data) to the popular (or not) vtk format. When time is specified the generated file will have a "_#" string in its name indicating its time order.

The program utilizes the PyVtk program and it should be used within the allowance range set by it.

=========== Execute the program ===========
To run the program simply use:

$ python vtk_converter.py output_file format_file file1 file2 ...

Here output_file is the name of the file for output without extensions. The format_file is a text file specifying the format for the convertion (explained later). The file list file1, file2, etc lists all the used files.


=========== The format file ===============
The format file specify how variables are formed from columns of data from file1, file2, etc. Each line of it should be of the form:
variable_name @ #1->#2, [#3->#4, #5->#6]
Here "variable_name" is the name of the variable (seen in e.g. paraview), "@" is a separater, #1 is the file index, and #2 is the column index. For example, if variable "ed" is stored in the 3rd column of "file1", this line should be written as:
ed @ 1 -> 3
The optional #3,4,5,6 play similar role as #1,2, and they are only needed if the variable being defined is a vector in which case each of them defines a component. #5,6 are not required and when skipped it is set to 0.
Four "variabl_name" string are special: x, y, z, and t. The first 3 defines the lattice that other data are attached to, and t is used to slice the data into time pieces. The omitted variables among x, y, and z will be automatically set to 0. If t is omitted, there will be only one file output, otherwise everytime t read from the data file is changing, a new file with a larger "_#" will be used. Typical output with time slicing looks like:
a_1.vtk, a_2.vtk, ...

Other special control parameters can also be specified in this file. For example, the line
skip @ 1->2, 2->5
means that the first 2 lines in file1 and first 5 lines in file2 will be skipped.


=========== Technical details =============
If the read (x,y,z) set form a "structured" grid, then a structuredGrid format will be used in the generated vtk file; otherwise a unstructuredGrid format will be used. See the convertToVtk function for details.

Data file must be purely numerical.

Enjoy.
