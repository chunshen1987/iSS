#! /usr/bin/env python
"""
    This script split the OSCAR.DAT file into multiple subfiles that 
    contain only small number of events. In this way, if there is a large
    number of events, we can run UrQMD using multiple CPU cores and also
    save a lot of re-run time if UrQMD crashes in some cases
"""

from sys import argv, exit
from os import makedirs, path

try:
    original_file = path.abspath(argv[1])
    number_of_split = int(argv[2])
except(IOError):
    print("Usage: split_events.py filename number_of_split")
    exit(1)

in_file = open(original_file, 'r')
out_file_list = []
for ifile in range(number_of_split):
    temp_file = open(original_file.split('.')[0] + "_%d.dat" % ifile, 'w')
    out_file_list.append(temp_file)

# header
for i in range(3):
    temp_line = in_file.readline()
    for ifile in range(number_of_split):
        out_file_list[ifile].write(temp_line)

file_idx = 0
while True:
    file_idx = file_idx % number_of_split
    event_line = in_file.readline()
    if event_line == '':
        break
    else:
        number_of_particles = int(event_line.split()[1])
    out_file_list[file_idx].write(event_line)
    for ipart in range(number_of_particles):
        temp_line = in_file.readline()
        out_file_list[file_idx].write(temp_line)
    file_idx += 1

in_file.close()
for ifile in range(number_of_split):
    out_file_list[ifile].close()
print("All events are split into %d files" % number_of_split)
