import pandas as pd
import numpy as np
import os
from argparse import ArgumentParser as AP


class DataFile(pd.DataFrame):
    def __init__(self, path, *args, **kwargs):
        super().__init__(self._harvest(path))
        self.path = os.path.realpath(path)
        self.name = path.split("/")[-1]
        self.xstep = self.iloc[1, 0] - self.iloc[0, 0]
        self.len = len(self.iloc[:, 0])

    def _harvest(self, path):
        """ Given the name of a file to read, harvest will read the data, return
        the column headers, and the (multicolumn) data as a numpy ndarray"""
        with open(path, 'r') as f:
            toprow = f.readline()
            try:
                float(toprow.split()[0])  # check to see if top row is headers
            except ValueError:  # top row is not headers
                head = 0
            else:
                head = None

        return (pd.read_csv(path, header=head))

    def _rescaleY(self, X, Y):
        """ Given the HHS from the dipole velocity or dipole length data, rescale it by
        omega^2 or omega^4 respectively to match the dipole acceleration"""
        if self.name.startswith("expec_v"):
            return Y*X*X
        elif self.name.startswith("expec_z"):
            return Y*X**4

    def _FFT(self, blackman=True, pad=8):
        """Apply a blackman window to all data columns, pad the data with
        leading zeros so that it is of length 2**pad, take fourier transform
        of the windowed input data, return DataFrame with transformed data"""
        df = pd.DataFrame()
        for col in self.columns[1:]:
            colname = "FT"+col
            ydat = self[col]
            if blackman:
                ydat = np.blackman(self.len)*ydat
            if pad:
                ydat = self._pad_with_zeros(ydat, factor=pad)
            Y = (np.fft.fft(ydat))[:len(ydat)//2]
            df[colname] = Y
        X = np.arange(len(Y))*(((np.pi)/len(Y))/self.xstep)
        df.insert(0, "Freq", X)

        return df

    def HHG(self):
        """Calculate the HHG spectrum"""
        df = self._FFT()
        X = df["Freq"]
        for col in df.columns[1:]:
            ydat = df[col]
            tmpdf = pd.DataFrame({col: ydat})
            df.update(tmpdf)
        return df

    def _pad_with_zeros(self, ydat, factor=8):
        """ pad ydat with zeros so that the length of ydat is a power of 2 and
        at least factor times the length of ydat on input."""

        pot = 2  # power of two
        while pot < factor * len(ydat):
            pot *= 2
        numzeros = pot-len(ydat)
        pre = (numzeros+1)//2
        post = numzeros//2
        return np.concatenate([np.zeros(pre), np.transpose(ydat),
                               np.zeros(post)])

# --------------------------------------------------------------------------------


def read_command_line():
    parser = AP()
    parser.add_argument('-f', '--dipoleFile',
                        help='path to csv file containing dipole '
                        'velocity output for various intentities',
                        default='expec_v.csv')
    parser.add_argument('-m', '--minimum_intensity', type=float,
                        help='change minimum intensity to use in '
                        'average, given as percentage of peak intensity, '
                        'rather than use all intensities in file')
    parser.add_argument('-di', '--delta_intensity', action='store_true',
                        help='option to calculate several averages '
                        'using various different values of delta I')
    parser.add_argument('-w', '--omega_0', action='store_true',
                        help='option to calculate several averages '
                        'using various different beam radii')
    return vars(parser.parse_args())

# --------------------------------------------------------------------------------


def get_intensities(df):
    "Gets the intensities used from the .csv file"
    intensities = np.array([float(x) for x in df.columns[1:]])
    peak_intensity = np.amax(intensities)
    return peak_intensity, intensities


def get_weights(I, I0, w0):
    "calculate the weights for gaussian laser beam in 2D configuration"
    return (np.pi*w0**2)/(2*I)


def phase(I, I0, w, w0, d=25):
    d *= 10**7  # convert from mm to nm
    w0 *= 10**6  # convert from mm to nm
    return (I/I0)**((1j*w0**2*w)/(4*d))


# --------------------------------------------------------------------------------
args = read_command_line()

df = DataFile(args['dipoleFile'])
peak_intensity, intensities = get_intensities(df)
dI = intensities[0]-intensities[1]

if args["minimum_intensity"]:
    # Get intensity corresponding to percentage of peak
    I_min = peak_intensity*args["minimum_intensity"]
    # Get intensity column corresponding to I_min at which to truncate
    I_trunc = intensities[np.abs(intensities-I_min).argmin()]
    df = df.truncate(before=str(I_trunc), after=None,
                     axis=1)  # truncate dataframe

vph = df.HHG()

# Write peak intensity spectra to file
# assumes peak intensity is first in csv file - generalise?
single = np.real(np.abs(vph.iloc[:, 1])**2)
single = df._rescaleY(vph["Freq"].values, single)

peak_hs = pd.DataFrame({"Freq": vph["Freq"].values*27.212, "0001_z": single})
peak_hs.to_csv("single_peak.csv", index=False)

# Loop over several different laser focal areas
areas = [0.01, 0.1, 1.5]
for area in areas:
    w0 = np.sqrt(area/np.pi)  # beam radius
    weights = get_weights(intensities, peak_intensity, w0)
    weights[-1] *= 0.5
    weights[0] *= 0.5
    vph = df.HHG()

    averaged_data_naive = 0
    averaged_data_full = 0

    for col, intensity, weight in zip(vph.columns[1:], intensities, weights):
        phase_factor = phase(intensity, peak_intensity, vph['Freq'].values, w0)
        averaged_data_naive += vph[col]*weight*dI
        averaged_data_full += vph[col]*phase_factor*weight*dI

    nc_avg = np.real(np.abs(averaged_data_naive)**2)
    nc_avg = df._rescaleY(vph["Freq"].values, nc_avg)

    # coherent phase weighted average
    fc_avg = np.real(
        np.abs((1j*vph["Freq"].values/(2*np.pi*25*10**7))*averaged_data_full)**2)
    fc_avg = df._rescaleY(vph["Freq"].values, fc_avg)

    # Write naive coherent average to file
    n_coh = pd.DataFrame({"Freq": vph["Freq"].values*27.212, "0001_z": nc_avg})
    n_coh.to_csv("naive_coherent_area_"+str(area)+"_step_" +
                 str(round(dI, 3))+".csv", index=False)

    # Write full coherent average to file
    f_coh = pd.DataFrame({"Freq": vph["Freq"].values*27.212, "0001_z": fc_avg})
    f_coh.to_csv("full_coherent_area_"+str(area)+"_step_" +
                 str(round(dI, 3))+".csv", index=False)


# Changing intensity step 2x bigger, 5x bigger, 10x bigger
for i in [2, 5, 10]:
    ndf = vph.iloc[:, 1::i]
    peak_intensity, intensities = get_intensities(df.iloc[:, 1::i])
    dI = intensities[0]-intensities[1]

    area = 0.1  # area of laser focus
    w0 = np.sqrt(area/np.pi)  # beam radius

    weights = get_weights(intensities, peak_intensity, w0)
    weights[-1] *= 0.5
    weights[0] *= 0.5

    averaged_data_naive = 0
    averaged_data_full = 0

    for col, intensity, weight in zip(ndf.columns, intensities, weights):
        phase_factor = phase(intensity, peak_intensity, vph['Freq'].values, w0)
        averaged_data_naive += ndf[col]*weight*dI
        averaged_data_full += ndf[col]*phase_factor*weight*dI

    nc_avg = np.real(np.abs(averaged_data_naive)**2)
    nc_avg = df._rescaleY(vph["Freq"].values, nc_avg)

    # coherent phase weighted average
    fc_avg = np.real(
        np.abs((1j*vph["Freq"].values/(2*np.pi*25*10**7))*averaged_data_full)**2)
    fc_avg = df._rescaleY(vph["Freq"].values, fc_avg)

    # Write naive coherent average to file
    n_coh = pd.DataFrame({"Freq": vph["Freq"].values*27.212, "0001_z": nc_avg})
    n_coh.to_csv("naive_coherent_area_"+str(area)+"_step_" +
                 str(round(dI, 3))+".csv", index=False)

    # Write full coherent average to file
    f_coh = pd.DataFrame({"Freq": vph["Freq"].values*27.212, "0001_z": fc_avg})
    f_coh.to_csv("full_coherent_area_"+str(area)+"_step_" +
                 str(round(dI, 3))+".csv", index=False)
