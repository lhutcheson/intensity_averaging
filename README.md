# Intensity Averaged high-harmonic spectra

This repository contains both the data and utilities for replicating the results found in the paper `Modelling intensity volume averaging in ab initio calculations of High Harmonic Generation`. 

The `bin/` directory contains the script used to calculate the intensity average, `gen_hhs_avg.py`. This reads in a csv file with each column containing the dipole expectation value at each laser intensity*.

This scripts output several different intensity averaged HHS using various different intensity increments and laser beam radiuses. The bin directory also contains scripts for plotting the figures used in the intensity averaging paper.

The data directory contains the RMT-generated dipole data used in the paper and the experimental data for comparison.

*For users of the RMT suite: a csv file in the correct format can be generated from rmt output files using the `bin/makedf.py` utility.

# Requirements
Python Packages:
	pandas
	numpy
	matplotlib
	argparse
	joblib
 
