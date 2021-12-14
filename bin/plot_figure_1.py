import numpy as np
import matplotlib.pyplot as plt

def moving_average(input_data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(input_data, window, 'same')

def smooth_data(wdata,nits=4):
    out_data = wdata
    for ii in np.arange(nits):
        out_data=moving_average(out_data,30)
    return out_data

def get_scale_factor(data_to_scale, energy_value):
   for i in range(len(data_to_scale[:,0])):
      if data_to_scale[i,0]>energy_value:
         normalisation_factor = data_to_scale[i,1]
         return(normalisation_factor)

exper = np.loadtxt('experimental_data.csv', skiprows=1, delimiter=',' )
data_1 = np.loadtxt('single_peak.csv', skiprows=1, delimiter=',')
data_2 = np.loadtxt('full_coherent_area_0.01_step_0.001.csv', skiprows=1, delimiter=',')
data_3 = np.loadtxt('full_coherent_area_0.1_step_0.001.csv', skiprows=1, delimiter=',')
data_4 = np.loadtxt('full_coherent_area_1.5_step_0.001.csv', skiprows=1, delimiter=',')

data_1[:,1] = smooth_data(data_1[:,1])
data_2[:,1] = smooth_data(data_2[:,1])
data_3[:,1] = smooth_data(data_3[:,1])
data_4[:,1] = smooth_data(data_4[:,1])

max_data_1 = get_scale_factor(data_1, 90)
max_data_2 = get_scale_factor(data_2, 90)
max_data_3 = get_scale_factor(data_3, 101)
max_data_4 = get_scale_factor(data_4, 90.8)
maxexper = get_scale_factor(exper, 90)

scalefactor= maxexper/maxexper
plt.semilogy(exper[:,0],scalefactor*exper[:,1], 'k', label="Experimental Data")

scale_factor = maxexper/max_data_1
plt.semilogy(data_1[:,0], scale_factor*data_1[:,1], 'r--', label='Single Intensity')
scale_factor = maxexper/max_data_2
plt.semilogy(data_2[:,0], scale_factor*data_2[:,1], label='$\omega_0=0.056$')
scale_factor = maxexper/max_data_3
plt.semilogy(data_3[:,0], scale_factor*data_3[:,1],label='$\omega_0=0.1784$')
scale_factor = maxexper/max_data_4
plt.semilogy(data_4[:,0], scale_factor*data_4[:,1],label='$\omega_0=0.6910$')

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')

plt.show()