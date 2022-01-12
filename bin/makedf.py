"""
makedf: wrangle RMT output data into the correct format for further postprocessing. 

This script reads the EField, expec_v_all and expec_z_all files from the
individual intensity directories generated with RMT, puts them together into
a dataframe, and outputs this to a .csv file. This csv file can then be read by
the intensity averaging utility to calculate various intensity averages.

The individual directories should be named after the intensity of the laser
pulse: intenity_X.XXX

The output file contains columns: "Time", "Intensity1", "Intensity2" ... where
"Intensity#" is the numerical value of the laser intensity in units of 10(14) Wcm(-2)
as per the required naming of the RMT directories (see previous paragraph)
"""
import pandas as pd
import glob
from argparse import ArgumentParser as AP


def read_command_line():
    parser = AP()
    parser.add_argument('files', nargs="+",
                        help="list of RMT output directories")
    parser.add_argument('-e', '--efield', action='store_true',
                        help='generate csv of all EField files')
    parser.add_argument('-v', '--expec_v', action='store_true',
                        help='generate csv of expec_v_all files')
    parser.add_argument('-z', '--expec_z', action='store_true',
                        help='generate csv of expec_z_all files')
    return vars(parser.parse_args())


def strip_num(filename):
    """strip numerical value of format X.XXX from the given filename"""
    import re
    pat = re.compile("[0-9]{1}.[0-9]{3}")
    m = pat.search(filename)
    return(float(m.group(0)))


def order_files(filelist):
    """put files in filelist in order of the numerical intensity values they
    represent"""
    from operator import itemgetter
    td_list = []
    for fname in filelist:
        td_list.append((fname, strip_num(fname)))
    return(sorted(td_list, key=itemgetter(1), reverse=True))


args = read_command_line()
filelist = args["files"]
filelist = order_files(filelist)

optdict = {"efield": "EField",
           "expec_v": "expec_v",
           "expec_z": "expec_z"}

for opt in optdict:
    if args[opt]:
        ndf = pd.DataFrame()
        for f, t in filelist: # f is filename, t is intensity value
            fname = glob.glob(f+f"/{optdict[opt]}*")[0]
            df = pd.read_csv(fname, delim_whitespace=True)
            ndf[str(t)] = df['0001_z']
        ndf.insert(0, "Time", df['Time'])
        ndf.drop(index=1)
        ndf.to_csv(f'{optdict[opt]}.csv', index=False)
