"""
This script reads the EField, expec_v_all and expec_z_all files from the individual 
intensity directories, and puts them together into a dataframe, which is output to 
a .csv file. This csv file can then be read by the intensity averaging utility to 
calculate various intensity averages.

The individual directories should be named after the intensity of the laser pulse: 
    intenity_X.XXX
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as sp
import glob
from argparse import ArgumentParser as AP

def read_command_line():
    parser = AP()
    parser.add_argument('files', nargs="+",help="list of directories containing the files for dataframe")
    parser.add_argument('-e','--efield', action='store_true', help='generate csv of all EField files')
    parser.add_argument('-v','--expec_v', action='store_true', help='generate csv of expec_v_all files')
    parser.add_argument('-z','--expec_z', action='store_true', help='generate csv of expec_z_all files')
    return vars(parser.parse_args())

def strip_num(filename):
    import re
    pat=re.compile("[0-9]{1}.[0-9]{3}")
    m=pat.search(filename)
    return(float(m.group(0)))

def order_files(filelist):
    from operator import itemgetter
    td_list=[]
    for fname in filelist:
        td_list.append((fname,strip_num(fname))) 
    return(sorted(td_list, key=itemgetter(1), reverse=True))

def grabData(fname):
    df = pd.read_csv(fname, delim_whitespace=True)
    return (df)

args = read_command_line()
filelist=args["files"]
filelist= order_files(filelist)

if args["efield"]:
    ndf = pd.DataFrame()
    for f,t in filelist:
        fname=glob.glob(f+"/EField*")[0]
        df = grabData(fname)
        ndf[str(t)] = df['0001_z']
    ndf.insert(0, "Time", df['Time'])
    ndf.drop(index=1)
    ndf.to_csv('EField.csv',index=False)

if args["expec_v"]:
    ndf = pd.DataFrame()
    for f,t in filelist:
        fname=glob.glob(f+"/expec_v_all*")[0]
        df = grabData(fname)
        ndf[str(t)] = df['0001_z']
    ndf.insert(0, "Time", df['Time'])
    ndf.drop(index=1)
    ndf.to_csv('expec_v.csv',index=False)

if args["expec_z"]:
    ndf = pd.DataFrame()
    for f,t in filelist:
        fname=glob.glob(f+"/expec_z_all*")[0]
        df = grabData(fname)
        ndf[str(t)] = df['0001_z']
    ndf.insert(0, "Time", df['Time'])
    ndf.drop(index=1)
    ndf.to_csv('expec_z.csv',index=False)


