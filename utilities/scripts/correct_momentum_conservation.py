#!/usr/bin/env python3

import sys

from numpy import *
from os import path

def parse_data(data_line):
    data = data_line.split()
    data = list(map(int, data[:2])) + list(map(float, data[2:]))
    return(data)

OSCAR_file_path = str(sys.argv[1])

OSCAR_file = open(OSCAR_file_path, 'r')
output_file = open('OSCAR_w_GMC.DAT', 'w')

line_count = 0
Nparticle = 0
event_header_line = 3
iev = 0
event_data = []
for temp_line in OSCAR_file:
    if line_count < 3:
        output_file.write(temp_line)

    if line_count == event_header_line:
        output_file.write(temp_line)
        Nparticle = int(temp_line.split()[1])
        event_header_line += Nparticle + 1
        iev += 1
        print("analysis event %d with %d particles" % (iev, Nparticle))
        event_data = []

    if line_count > (event_header_line - Nparticle - 1):
        data = parse_data(temp_line)
        event_data.append(data)

    if line_count > 3 and line_count == event_header_line - 1:
        event_data = array(event_data)
        CM_Px = sum(event_data[:, 2])
        CM_Py = sum(event_data[:, 3])
        CM_Pz = sum(event_data[:, 4])
        total_E = sum(event_data[:, 5])
        print("total energy = %g GeV" % total_E)
        print("correction per particle: delta_px = %g GeV, "
              "delta_py = %g GeV, delta_pz = %g GeV"
              % (CM_Px/Nparticle, CM_Py/Nparticle, CM_Pz/Nparticle))
        event_data[:, 2] -= CM_Px/Nparticle
        event_data[:, 3] -= CM_Py/Nparticle
        event_data[:, 4] -= CM_Pz/Nparticle
        event_data[:, 5] = sqrt(event_data[:, 2]**2. + event_data[:, 3]**2.
                                + event_data[:, 4]**2. + event_data[:, 6]**2.)
        for iline in range(Nparticle):
            output_file.write(
                "%d  %d  %.16e  %.16e  %.16e  %.16e  %.16e  %.16e  %.16e  %.16e  %.16e\n"
                % (event_data[iline, 0], event_data[iline, 1],
                   event_data[iline, 2], event_data[iline, 3],
                   event_data[iline, 4], event_data[iline, 5],
                   event_data[iline, 6], event_data[iline, 7],
                   event_data[iline, 8], event_data[iline, 9],
                   event_data[iline, 10]))

    line_count += 1

OSCAR_file.close()
output_file.close()
