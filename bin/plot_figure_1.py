"""
plot figure 1 of paper
"""

import numpy as np
import matplotlib.pyplot as plt

def moving_average(input_data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(input_data, window, 'same')

def smooth_data(wdata,nits=8):
    out_data = wdata
    for ii in np.arange(nits):
        out_data=moving_average(out_data,80)
    return out_data

def get_scale_factor(data_to_scale, energy_value=101, col=-1):
   for i in range(len(data_to_scale[:,0])):
      if data_to_scale[i,0]>energy_value:
         normalisation_factor = data_to_scale[i,col]
         return(normalisation_factor)

def get_scale_factor(data_to_scale, energy_value=101, col=-1):
   for i in range(len(data_to_scale[:,0])):
      if data_to_scale[i,0]>energy_value:
         normalisation_factor = data_to_scale[i,col]
         return(max_exper/normalisation_factor)

exper = np.loadtxt('experimental_data.csv', skiprows=1, delimiter=',' )
max_exper = get_exp_scale_factor(exper)
plt.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data")

files=['single_peak.csv', 'naive_coherent_area_0.1_step_1.csv', 'intermediate_coherent_area_0.1_step_1.csv', 'full_coherent_area_0.1_step_1.csv',]
labels=['Single Intensity, $I_{max} = 1.9$x$10^{14}$ W cm$^{-2}$', 'Naive Model', 'Intermediate Model', 'Full Model']

for file, label in zip(files,labels):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,-1] = smooth_data(data[:,-1])
   scale_factor = max_exper/get_scale_factor(data)
   plt.semilogy(data[:,0], scale_factor*data[:,-1], linewidth='1.5', label=label)

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')
plt.show()