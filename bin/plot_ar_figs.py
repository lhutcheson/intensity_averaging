"""
Script for plotting the Argon results in:
   'Modelling intensity volume averaging in ab initio 
    calculations of High Harmonic Generation'
The data required for plotting may be produced using the
`gen_hhs_avg.py` script found on the repository, and two
different ranges of intensities are used:
   I_min = 0.7 I_max, and
   I_min = 0.95 I_max.

"""

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


def smooth_data(wdata,window_size=150,nits=30):
    out_data = wdata
    for ii in np.arange(nits):
        out_data=moving_average(out_data,window_size)
    return out_data

fig, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.set_xlim([20,102])
ax1.set_ylim([10e-1, 10e6])
ax2.set_ylim([10e-1, 10e6])

exper = np.loadtxt('experimental_data.csv', skiprows=1)
single = np.loadtxt('all_harm.csv', skiprows=1, delimiter=',')

max_exper = get_scale_factor(exper, energy_value=75)
scale_factor = max_exper/get_scale_factor(single, energy_value=75, col=1)
single[:,1] *= scale_factor
single[:,1] = smooth_data(single[:,1], window_size=60, nits=10)

ax2.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data: $3.5x 10^{14} W cm^{-2}$")
ax1.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5')

ax2.semilogy(single[:,0],single[:,1], 'r--', linewidth='1.5', label="Single Intensity: $4.1x 10^{14} W cm^{-2}$")
ax1.semilogy(single[:,0],single[:,1], 'r--', linewidth='1.5')

# ==================================================================
# First set of data: Imin=0.7Imax
# ==================================================================

files=['I_min0.7/naive_step_1.csv',
         'I_min0.7/full_area_0.05_step_1.csv'
         ]


labels=['Naive Method',
         'Full Method'
         ]

points = [75, 75]

for file, label, point in zip(files,labels, points):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,1] = smooth_data(data[:,1])
   max_exper = get_scale_factor(exper, energy_value=point)
   scale_factor = max_exper/get_scale_factor(data, energy_value=point, col=1)
   ax1.semilogy(data[:,0], scale_factor*data[:,1], linewidth='1.5', label=label)

# ==================================================================
# Second set of data: Imin=0.9Imax
# ==================================================================

files=[
         'I_min0.95/naive_step_1.csv',
         'I_min0.95/full_area_0.01_step_1.csv'
         ]


labels=[
         'Naive Method',
         'Full Method'
         ]

points = [75, 82.9]

for file, label, point in zip(files,labels, points):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,1] = smooth_data(data[:,1])
   max_exper = get_scale_factor(exper, energy_value=point)
   scale_factor = max_exper/get_scale_factor(data, energy_value=point, col=1)
   ax2.semilogy(data[:,0], scale_factor*data[:,1],linewidth='1.5', label=label)

ax2.legend(fontsize="x-small")
fig.supxlabel('Energy (eV)')
fig.supylabel('Harmonic Intensity (arb.)')
# ax1.set_title('$I_{min} = 0.70 I_{max}$')
# ax2.set_title('$I_{min} = 0.90 I_{max}$')
plt.tight_layout()
plt.show()

