import numpy as np
import matplotlib.pyplot as plt

def moving_average(input_data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(input_data, window, 'same')

def smooth_data(wdata,nits=10):
    out_data = wdata
    for ii in np.arange(nits):
        out_data=moving_average(out_data,80)
    return out_data

def get_scale_factor(data_to_scale, energy_value=100, col=-1):
   for i in range(len(data_to_scale[:,0])):
      if data_to_scale[i,0]>energy_value:
         normalisation_factor = data_to_scale[i,col]
         return(normalisation_factor)


exper = np.loadtxt('experimental_data.csv', skiprows=1, delimiter=',' )
max_exper = get_scale_factor(exper)
plt.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data")

files=['single_peak.csv', 'full_coherent_area_0.01_step_1.csv', 'full_coherent_area_0.1_step_1.csv', 'full_coherent_area_1.5_step_1.csv',]
labels=['Single Intensity, $I_{max} = 1.9$x$10^{14}$ W cm$^{-2}$', '$\omega_0=0.056$mm', '$\omega_0=0.1784$mm', '$\omega_0=0.6910$mm']

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

# exper = np.loadtxt('experimental_data.csv', skiprows=1, delimiter=',' )
# data_1 = np.loadtxt('single_peak.csv', skiprows=1, delimiter=',')
# data_2 = np.loadtxt('full_coherent_area_0.01_step_1.csv', skiprows=1, delimiter=',')
# data_3 = np.loadtxt('full_coherent_area_0.1_step_1.csv', skiprows=1, delimiter=',')
# data_4 = np.loadtxt('full_coherent_area_1.5_step_1.csv', skiprows=1, delimiter=',')

# data_1[:,1] = smooth_data(data_1[:,1])
# data_2[:,-1] = smooth_data(data_2[:,-1])
# data_3[:,-1] = smooth_data(data_3[:,-1])
# data_4[:,-1] = smooth_data(data_4[:,-1])

# maxexper = get_scale_factor(exper, 90)
# max_data_1 = get_scale_factor(data_1, 90)
# max_data_2 = get_scale_factor(data_2, 100)
# max_data_3 = get_scale_factor(data_3, 100)
# max_data_4 = get_scale_factor(data_4, 100)

# scalefactor= maxexper/maxexper
# plt.semilogy(exper[:,0],scalefactor*exper[:,1], 'k', linewidth='2.5', label="Experimental Data")

# scale_factor = maxexper/max_data_1
# plt.semilogy(data_1[:,0], scale_factor*data_1[:,1], 'r', linewidth='1.5', label='Single Intensity')

# maxexper = get_scale_factor(exper, 100)
# scale_factor = maxexper/max_data_2
# plt.semilogy(data_2[:,0], scale_factor*data_2[:,-1], linestyle=':', linewidth='1.5', label='$\omega_0=0.056$mm')
# scale_factor = maxexper/max_data_3
# plt.semilogy(data_3[:,0], scale_factor*data_3[:,-1], linestyle='dashdot', linewidth='1.5', label='$\omega_0=0.1784$mm')
# scale_factor = maxexper/max_data_4
# plt.semilogy(data_4[:,0], scale_factor*data_4[:,-1], linestyle='dashed', linewidth='1.5',label='$\omega_0=0.6910$mm')

# plt.xlabel('Photon Energy (eV)')
# plt.ylabel('Harmonic Intensity (arb. units)')
# plt.xlim([0,160])
# plt.ylim([10e-6, 10e-2])
# plt.legend(fontsize='small')

# plt.show()