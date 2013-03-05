#!/usr/bin/env python

# system-wide
from sys import argv
from os import path
from os import getcwd
from time import strftime

# my libraries
from runR import run
from fileR import makeDir

def record_spectra(argv):
  """ Run the commands specified by argv[2:] with freeze-out data stored
  under argv[1]/results. """

  if len(argv) > 1:
    # generate necessary strings
    unique_signature = "-".join(argv[2:])+"--"+strftime("%Y-%m-%d@%H:%M:%S")
    logfile = unique_signature+".log"
    exe_line = " ".join(argv[2:])+" | tee "+logfile

    # clean results folder for spectra calculation
    run("bash ./prepare-spectra.sh"+" "+argv[1])

    # execute hydro code (presumably)
    print "Executing " + exe_line + "..."
    run(exe_line)

    # copy out generated data files
    run("bash ./copy-data-spectra.sh"+" "+"."+" "+argv[1])

  else:
    print "Usage: record-spectra dir-with-hydro-results command (with parameters)"

if __name__ == "__main__":
  record_spectra(argv)
