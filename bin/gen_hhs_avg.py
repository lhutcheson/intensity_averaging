import pandas as pd
import numpy as np
import os
from argparse import ArgumentParser as AP
import matplotlib.pyplot as plt
from joblib import Parallel, delayed


class DipoleFile(pd.DataFrame):
    """ Class for handling data file of time-dependent expectation value of the
    dipole (d(t)). In particular, this class handles harvesting the data from
    file or existing pandas dataFrame, santitising the data, and fourier
    transforming it to give the harmonic spectrum. Parameter "form" corresponds
    either to "v" velocity form, or "l" length form of the dipole.
    """

    def __init__(self, path=None, dF=None, form="v", *args, **kwargs):
        if dF:
            super().__init__(dF)
        elif path:
            super().__init__(self._harvest(path))
        self.xstep = self.iloc[1, 0] - self.iloc[0, 0]
        self.len = len(self.iloc[:, 0])
        self._intensities()  # sets peak_intensity
        if form == "v":
            self.scalefactor = 2
        else:
            self.scalefactor = 4

    def _harvest(self, path):
        """ Given the name of a file to read, harvest will read the data, return
        the column headers, and the (multicolumn) data as a pandas DataFrame"""
        with open(path, 'r') as f:
            toprow = f.readline()
            try:
                float(toprow.split()[0])  # check to see if top row is headers
            except ValueError:  # top row is not headers
                head = 0
            else:
                head = None

        return (pd.read_csv(path, header=head))

    def _intensities(self):
        """Gets the intensities used from the .csv file, units 10^14 Wcm^-2"""
        intens = [float(x) for x in self.columns[1:]]
        intens.sort(reverse=True)  # ensures descending order
        intens = np.array(intens)
        self.peak_intensity = np.amax(intens)
        return intens

    def _FFT(self, blackman=True, pad=8, collim=None):
        """Apply a blackman window to all data columns, pad the data with
        leading zeros so that it is of length 2**pad, take fourier transform
        of the windowed input data, return DataFrame with transformed data"""
        df = pd.DataFrame()
        if collim is None:
            collist = self.columns[1:]
        else:
            collist = self.columns[1:collim+1]
        for col in collist:
            ydat = self[col]
            if blackman:
                ydat = np.blackman(self.len)*ydat
            if pad:
                ydat = self._pad_with_zeros(ydat, factor=pad)
            Y = (np.fft.fft(ydat))[:len(ydat)//2]
            df[col] = Y
        X = np.arange(len(Y))*(((np.pi)/len(Y))/self.xstep)
        df.insert(0, "Freq", X)

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

    def peakHHG(self):
        df = self._FFT(collim=1)
        col = df.columns[1]
        amplitude = np.real(np.abs(df[col].values)**2)
        amplitude = amplitude * df["Freq"]**self.scalefactor
        df["Freq"] *= 27.212  # convert frequency to eV
        df[col] = amplitude
        newdf = df[["Freq", col]]
        return newdf

    def _weights(self, I, I0, w0):
        """calculate the weights for gaussian laser beam in 2D configuration"""
        weights = (np.pi*w0**2)/(2*I)
        weights[0] *= 0.5
        weights[-1] *= 0.5
        return weights

    def _phase(self, I, I0, w, w0, d=25):
        d *= 10**7  # convert from cm to nm
        w0 *= 10**6  # convert from mm to nm
        return (I/I0)**((1j*w0**2*w)/(4*d))

    def _detector(self, R_detect, I, I0, w0, d=25, w=np.array([4.55949])):
        """Account for the different detection points, R_detect, of the detector"""
        import scipy.special as sp

        R_detect *= 10**6  # convert from mm to nm
        w0 *= 10**6  # convert from mm to nm
        d *= 10**7  # convert from cm to nm
        L = 45.5949/w # convert omega in a.u to wavelength in nm

        r_sq = -0.5*w0**2*np.log(I/I0)

        a = (-2*np.pi*np.sqrt(r_sq)*R_detect)/(d*L)
        a = round(max(abs(a)))
        pre_factor = np.exp((np.pi*1j*(R_detect**2 + r_sq))/(d*L))

        term =0
        for l in range(0, 2*a, 2):
            term += (2*l+1)*(1j)**l*sp.spherical_jn(l,a)*((np.pi/4)**l)*(sp.binom(l,l/2))**2

        return 2*pre_factor*term

    def gen_radius_contribution(self, radius, df, stride, intensities, weights, I0, w, w0, d=25):
        """calculate the contribution for the different detection points, R_detect, of the detector"""
        dI = intensities[1]-intensities[0]
        int_summed = 0
        for col, intensity, weight in zip(
            df.columns[1::stride], intensities, weights):
            phase_factor = self._phase(
                intensity, I0, w, w0, d)
            detector_factor = self._detector(
                radius, intensity, I0, w0, d, w)
            int_summed += detector_factor * phase_factor * weight * dI * df[col]
        
        pre_factor = (1j*w)/(2*np.pi*d*10**7)
        int_summed *= pre_factor
        
        return int_summed


    def intensityAveragedHHG(self, no_cores, focus=0.1, d=25, R_detect=None,
                             phase=True, stride=1, I_min=None):
        self._select(I_min)
        df = self._FFT()
        w = df["Freq"]
        w0 = np.sqrt(focus/np.pi)  # beam radius
        intensities = self._intensities()[::stride]
        dI = intensities[0] - intensities[1]
        weights = self._weights(intensities, self.peak_intensity, w0)
        peak = self.peak_intensity
        if R_detect:
            rad_summed = 0

            radiuses = np.arange(0.0001,R_detect+1e-06,1e-06)
            radius_contributions = Parallel(n_jobs=no_cores)(delayed(self.gen_radius_contribution)(radius, df, stride, intensities, weights, peak, w, w0) for radius in radiuses)
            
            coh_avg = np.real(np.abs(sum(radius_contributions))**2)*w**self.scalefactor
            incoh_avg = sum(np.real(np.abs(cont)**2)*w**self.scalefactor for cont in radius_contributions)

            rad_avg = pd.DataFrame({"Freq": w*27.212, "incoh_avg": incoh_avg, "coh_avg": coh_avg})

            return rad_avg

        else:
            summed = 0
            for col, intensity, weight in zip(
                df.columns[1::stride], intensities, weights):
                if phase:
                    phase_factor = self._phase(
                        intensity, self.peak_intensity, w, w0, d)
                else:
                    phase_factor = 1.0
                summed += phase_factor * weight * dI * df[col]
            
            if phase:
                pre_factor = (1j*w)/(2*np.pi*d*10**7)
            else:
                pre_factor = 1.0

            summed *= pre_factor

            avg = np.real(np.abs(summed)**2)*w**self.scalefactor
            newdf = pd.DataFrame({"Freq": w*27.212, "0001_z": avg})
        
            return newdf            

    def _select(self, minimum_intensity):
        """If a minimum intensity threshold is set, drop all columns in the
        DipoleFile which correspond to intensities lower than this threshold"""
        if minimum_intensity:
            # Get intensity corresponding to percentage of peak
            I_min = self.peak_intensity*args["minimum_intensity"]
            print('Imin: ', I_min)
            # truncate dataframe
            for col in self.columns[1:]:
                columnIntensity = float(col)
                print('Intensity: ', columnIntensity)
                if columnIntensity < I_min:
                    self.drop(col, axis=1, inplace=True)


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
    parser.add_argument('-p', '--parallel', type=int, default=1,
                        help='option to run the code in parallel'
                        'and select the number of cores used')
    return vars(parser.parse_args())

# --------------------------------------------------------------------------------

if __name__ == "__main__":
    args = read_command_line()
    I_min = args["minimum_intensity"]
    no_cores = args['parallel']

    df = DipoleFile(args['dipoleFile'])


# Write peak intensity spectra to file
    # single = df.peakHHG()
    # single.to_csv("single_peak.csv", index=False)

# Loop over several different laser focal areas
    for area in [0.1]:#[0.01, 0.1, 1.5]:
        # Write naive coherent average to file
        # n_coh = df.intensityAveragedHHG(focus=area, d=25, phase=False, I_min=I_min)
        # n_coh.to_csv(f"naive_coherent_area_{area}_step_1.csv",
        #              index=False)

        # Write intermediate coherent average to file
        # f_coh = df.intensityAveragedHHG(focus=area, d=25, phase=True, I_min=I_min)
        # f_coh.to_csv(f"intermediate_coherent_area_{area}_step_1.csv",
        #              index=False)

        # Write full coherent average accounting for detector to file
        f_coh_det = df.intensityAveragedHHG(no_cores=no_cores, focus=area, R_detect=0.0005, d=25, phase=True, I_min=I_min)
        f_coh_det.to_csv(f"full_coherent_area_{area}_step_1.csv",
                     index=False)

# Changing intensity step 2x bigger, 5x bigger, 10x bigger
    # area=0.1
    # for stride in [2, 5, 10]:
    #     df = DipoleFile(args['dipoleFile'])

    #     n_coh = df.intensityAveragedHHG(focus=area, d=25, phase=False, I_min=I_min,
    #                                     stride=stride)
    #     n_coh.to_csv(f"naive_coherent_area_0.1_step_{stride}.csv",
    #                  index=False)

        # # Write intermediate coherent average to file
        # f_coh = df.intensityAveragedHHG(focus=area, d=25, phase=True, I_min=I_min,
        #                                 stride=stride)
        # f_coh.to_csv(f"intermediate_coherent_area_0.1_step_{stride}.csv",
        #              index=False)

    #     # Write full coherent average accounting for detector to file
        # f_coh_det = df.intensityAveragedHHG(focus=area, d=25, R_detect=0.0005, 
        #                                     phase=True, I_min=I_min, stride=stride)
        # f_coh_det.to_csv(f"full_coherent_area_0.1_step_{stride}.csv",
        #              index=False)

