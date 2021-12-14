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
data_2 = np.loadtxt('naive_coherent_area_0.1_step_0.001.csv', skiprows=1, delimiter=',')
data_3 = np.loadtxt('full_coherent_area_0.1_step_0.001.csv', skiprows=1, delimiter=',')

smooth_1 = np.loadtxt('single_peak.csv', skiprows=1, delimiter=',')
smooth_2 = np.loadtxt('naive_coherent_area_0.1_step_0.001.csv', skiprows=1, delimiter=',')
smooth_3 = np.loadtxt('full_coherent_area_0.1_step_0.001.csv', skiprows=1, delimiter=',')
smooth_1[:,1] = smooth_data(smooth_1[:,1])
smooth_2[:,1] = smooth_data(smooth_2[:,1])
smooth_3[:,1] = smooth_data(smooth_3[:,1])

max_data_1 = get_scale_factor(data_1, 64.68)
max_data_2 = get_scale_factor(data_2, 64.68)
max_data_3 = get_scale_factor(data_3, 65.4)
max_smooth_1 = get_scale_factor(smooth_1, 64.68)
max_smooth_2 = get_scale_factor(smooth_2, 64.68)
max_smooth_3 = get_scale_factor(smooth_3, 65.4)
maxexper = get_scale_factor(exper, 64.68)

scalefactor= maxexper/maxexper
plt.semilogy(exper[:,0],scalefactor*exper[:,1], 'k', label="Experimental Data")

scale_factor = maxexper/max_data_1
plt.semilogy(data_1[:,0], scale_factor*data_1[:,1], '#1f77b4', alpha=0.5)
scale_factor = maxexper/max_data_2
plt.semilogy(data_2[:,0], scale_factor*data_2[:,1], '#ff7f0e', alpha=0.5)
scale_factor = maxexper/max_data_3
plt.semilogy(data_3[:,0], scale_factor*data_3[:,1], '#2ca02c', alpha=0.5)

scale_factor = maxexper/max_smooth_1
plt.semilogy(smooth_1[:,0], scale_factor*smooth_1[:,1], '#1f77b4', label='Single Intensity, $I_{max} = 6.58$x$10^{14}$ W cm$^{-2}$)')
scale_factor = maxexper/max_smooth_2
plt.semilogy(smooth_2[:,0], scale_factor*smooth_2[:,1], '#ff7f0e', label='Naive Model')
scale_factor = maxexper/max_smooth_3
plt.semilogy(smooth_3[:,0], scale_factor*smooth_3[:,1], '#2ca02c',label='Full Model')

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-1])
plt.legend(fontsize='small')

plt.show()