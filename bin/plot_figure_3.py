"""
Utility to plot Figure 3 of paper
"""
import numpy as np
import matplotlib.pyplot as plt

def moving_average(input_data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(input_data, window, 'same')

def smooth_data(wdata,nits=4):
    out_data = wdata
    for ii in np.arange(nits):
        out_data=moving_average(out_data,60)
    return out_data

def get_scale_factor(data_to_scale, energy_value):
   for i in range(len(data_to_scale[:,0])):
      if data_to_scale[i,0]>energy_value:
         normalisation_factor = data_to_scale[i,-1]
         return(normalisation_factor)

exper = np.loadtxt('experimental_data.csv', skiprows=1, delimiter=',')
max_exper = get_scale_factor(exper, 90)

naive_list = ['naive_coherent_area_0.1_step_1.csv', 'naive_coherent_area_0.1_step_2.csv', 'naive_coherent_area_0.1_step_5.csv', 'naive_coherent_area_0.1_step_10.csv']
intermediate_list = ['intermediate_coherent_area_0.1_step_1.csv', 'intermediate_coherent_area_0.1_step_2.csv', 'intermediate_coherent_area_0.1_step_5.csv', 'intermediate_coherent_area_0.1_step_10.csv']
full_list = ['full_coherent_detector_area_0.1_step_1.csv', 'full_coherent_detector_area_0.1_step_2.csv', 'full_coherent_detector_area_0.1_step_5.csv', 'full_coherent_detector_area_0.1_step_10.csv']
step_list = ['1', '2', '5', '10']

plt.figure(1)
plt.semilogy(exper[:,0], exper[:,1], 'k', label="Experimental Data")

for i in step_list:
   smooth = np.loadtxt(f'naive_coherent_area_0.1_step_{i}.csv', skiprows=1, delimiter=',')
   smooth[:,1] = smooth_data(smooth[:,1])
   max_smooth = get_scale_factor(smooth, 90)
   scale_factor = max_exper/max_smooth
   plt.semilogy(smooth[:,0], scale_factor*smooth[:,1], label=f'$\Delta I = ${i}$x10^{11}$')
plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')

plt.figure(2)
plt.semilogy(exper[:,0], exper[:,1], 'k', label="Experimental Data")
max_exper = get_scale_factor(exper, 101)

for i in step_list:
   smooth = np.loadtxt(f'intermediate_coherent_area_0.1_step_{i}.csv', skiprows=1, delimiter=',')
   smooth[:,1] = smooth_data(smooth[:,1])
   max_smooth = get_scale_factor(smooth, 101)

   scale_factor = max_exper/max_smooth
   plt.semilogy(smooth[:,0], scale_factor*smooth[:,1], label=f'$\Delta I = ${i}$x10^{11}$')

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')

plt.figure(3)


plt.semilogy(exper[:,0], exper[:,1], 'k', label="Experimental Data")
max_exper = get_scale_factor(exper, 101)

for i in step_list:
   smooth = np.loadtxt(f'full_coherent_area_0.1_step_{i}.csv', skiprows=1, delimiter=',')
   smooth[:,-1] = smooth_data(smooth[:,-1])
   max_smooth = get_scale_factor(smooth, 101)

   scale_factor = max_exper/max_smooth
   plt.semilogy(smooth[:,0], scale_factor*smooth[:,-1], label=f'$\Delta I = ${i}$x10^{11}$')

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')

plt.show()

