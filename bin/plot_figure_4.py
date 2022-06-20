import numpy as np
import matplotlib.pyplot as plt

def get_scale_factor(data_to_scale, energy_value, col=-1):
   for i in range(len(data_to_scale[:,0])):
      if data_to_scale[i,0]>energy_value:
         normalisation_factor = data_to_scale[i,col]
         return(normalisation_factor)
         
def moving_average(input_data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(input_data, window, 'same')

def smooth_data(wdata,nits=4):
    out_data = wdata
    for ii in np.arange(nits):
        out_data=moving_average(out_data,60)
    return out_data

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.set_xlim([20,102])
ax1.set_ylim([10e-6, 10e5])
ax2.set_ylim([10e-6, 10e5])
ax3.set_ylim([10e-6, 10e5])

# ==================================================================
# First set of data: Exp. intensity = 2.5e14
# ==================================================================

exper = np.loadtxt('./int_2.5/experimental_data.csv', skiprows=1)
energy_value = 65.12
max_exper = get_scale_factor(exper, energy_value)
ax1.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data: $2.5x 10^{14} W/cm^2$")

files=['./int_2.5/single_peak.csv', './int_2.5/naive_coherent_area_0.1_step_1.csv', './int_2.5/full_coherent_area_0.1_step_1.csv']
labels=["Single Intensity: $3.1x 10^{14} W/cm^2$", "Naive Model", "Full Model"]

for file, label in zip(files,labels):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,-1] = smooth_data(data[:,-1])
   scale_factor = max_exper/get_scale_factor(data, energy_value)
   ax1.semilogy(data[:,0], scale_factor*data[:,-1], linewidth='1.5', label=label)

ax1.legend(loc='lower left', fontsize="x-small")

# ==================================================================
# Second set of data: Exp. intensity = 2.9e14
# ==================================================================

exper = np.loadtxt('./int_2.9/experimental_data.csv', skiprows=1)
energy_value = 76.2
max_exper = get_scale_factor(exper, energy_value)
ax2.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data: $2.9x 10^{14} W/cm^2$")

files=['./int_2.9/single_peak.csv', './int_2.9/naive_coherent_area_0.1_step_1.csv', './int_2.9/full_coherent_area_0.1_step_1.csv']
labels=["Single Intensity: $3.5x 10^{14} W/cm^2$", "Naive Model", "Full Model"]

for file, label in zip(files,labels):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,-1] = smooth_data(data[:,-1])
   scale_factor = max_exper/get_scale_factor(data, energy_value)
   ax2.semilogy(data[:,0], scale_factor*data[:,-1], linewidth='1.5', label=label)

ax2.legend(loc='lower left', fontsize="x-small")

# ==================================================================
# Third set of data: Exp. intensity = 3.5e14
# ==================================================================

exper = np.loadtxt('./int_3.5/experimental_data.csv', skiprows=1)
energy_value = 85.9
max_exper = get_scale_factor(exper, energy_value)
ax3.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data: $3.5x 10^{14} W/cm^2$")

files=['./int_3.5/single_peak.csv', './int_3.5/naive_coherent_area_0.1_step_1.csv', './int_3.5/full_coherent_area_0.1_step_1.csv']
labels=["Single Intensity: $4.1x 10^{14} W/cm^2$", "Naive Model", "Full Model"]

for file, label in zip(files,labels):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,-1] = smooth_data(data[:,-1])
   scale_factor = max_exper/get_scale_factor(data, energy_value)
   ax3.semilogy(data[:,0], scale_factor*data[:,-1], linewidth='1.5', label=label)

ax3.legend(loc='lower left', fontsize="x-small")

plt.xlabel('Energy (eV)')
ax2.set_ylabel('Harmonic Intensity (arb.)')
plt.tight_layout()
plt.show()

