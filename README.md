# Intensity Averaged high-harmonic spectra

This repository contains utilities for computing intensity averaged high-harmonic spectra from dipoles calculated at single intensities.

The `bin/` directory contains the script used to calculate the intensity average, `gen_hhs_avg.py`. This reads in a csv file with each column containing the dipole expectation value at each laser intensity*.

This scripts output several different intensity averaged HHS using various different intensity increments and laser beam radii. The bin directory also contains scripts for plotting the figures used in the intensity averaging paper.

The data directory contains the data used in the paper.

*For users of the RMT suite: a csv file in the correct format can be generated from rmt output files using the `bin/makedf.py` utility. 
