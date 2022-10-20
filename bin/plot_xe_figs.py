"""
Script for plotting the Xenon results in:
   'Modelling intensity volume averaging in ab initio 
    calculations of High Harmonic Generation'
The data required for plotting may be produced using the
`gen_hhs_avg.py` script found on the repository, and two
different ranges of intensities are used:
   I_min = 0.6 I_max, and
   I_min = 0.95 I_max.

"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd

def moving_average(input_data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(input_data, window, 'same')

def smooth_data(wdata,window_size=100,nits=20):
    out_data = wdata
    for ii in np.arange(nits):
        out_data=moving_average(out_data,window_size)
    return out_data

def get_scale_factor(data_to_scale, energy_value=101, col=-1):
   for i in range(len(data_to_scale[:,0])):
      if data_to_scale[i,0]>energy_value:
         normalisation_factor = data_to_scale[i,col]
         return(normalisation_factor)

def get_scale_factor_alt(energy, data_to_scale, energy_value=101):
   for i in range(len(data_to_scale)):
      if energy[i]>energy_value:
         normalisation_factor = data_to_scale[i]
         return(normalisation_factor)


exper = np.loadtxt('experimental_data.csv', skiprows=1, delimiter=',' )
max_exper = get_scale_factor(exper)


# =====================================================


plt.figure(1)
# plt.title('Figure 1 from paper')
plt.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data")

files=[
         'Imin_95/single_peak.csv',
         'Imin_95/naive_step_1.csv',
         'Imin_95/full_area_0.05_step_1.csv',
         'Imin_60/naive_step_1.csv',
         'Imin_60/full_area_0.01_step_1.csv'
         ]


labels=[
         'Peak Intensity',
         'Naive, I$_{min}$ = 95% I$_{max}$',
         'Full, I$_{min}$ = 95% I$_{max}$',
         'Naive, I$_{min}$ = 60% I$_{max}$',
         'Full, I$_{min}$ = 60% I$_{max}$'
         ]

norm_points=[101, 110, 72, 110, 72]
cols = [-1,-1,-2,-1,-2]
styles = ['r--','-','-','-','-']
window_size = [100,100,100,300,300]
for file, label, point, col, style, window in zip(files,labels, norm_points, cols, styles, window_size):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,col] = smooth_data(data[:,col], window)
   max_exper = get_scale_factor(exper, energy_value=point)
   scale_factor = max_exper/get_scale_factor(data, energy_value=point, col=col)
   plt.semilogy(data[:,0], scale_factor*data[:,col], style, linewidth='1.5', label=label)

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')
plt.tight_layout()


# ============================================================

plt.figure(2)
# plt.title('Single intensity spectra')

file = pd.read_csv('all_harm.csv')

cols = [1.9, 1.711, 1.519, 1.33, 1.141, 0.949, 0.76, 0.568]
labels=[ '$I_{max}$',
         '$0.90 I_{max}$',
         '$0.80 I_{max}$',
         '$0.70 I_{max}$',
         '$0.60 I_{max}$',
         '$0.50 I_{max}$',
         '$0.40 I_{max}$',
         '$0.30 I_{max}$'
         ]

for col, label in zip(cols,labels):
   plt.semilogy(file['Freq'], smooth_data(file[str(col)]), linewidth='1.5', label=label)

max_exper = get_scale_factor(exper, energy_value=101)
max_data = get_scale_factor_alt(file['Freq'], file['1.9'], energy_value=101)

plt.semilogy(exper[:,0],(max_data/max_exper)*exper[:,1], 'k', linewidth='1.5', label="Experimental Data")

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')
plt.tight_layout()

# =====================================================


plt.figure(3)
# plt.title('Figure 3 from paper')
plt.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data")

files=[
         'Imin_95/naive_step_1.csv',
         'Imin_95/naive_step_2.csv',
         'Imin_95/naive_step_5.csv',
         'Imin_95/naive_step_10.csv',
         'Imin_95/naive_step_12.csv'
         ]


labels=[
         '$\Delta$I = 1x10$^11$ Wcm$^{-2}$',
         '$\Delta$I = 2x10$^11$ Wcm$^{-2}$',
         '$\Delta$I = 5x10$^11$ Wcm$^{-2}$',
         '$\Delta$I = 10x10$^11$ Wcm$^{-2}$',
         '$\Delta$I = 12x10$^11$ Wcm$^{-2}$',
         ]

norm_points=[101, 101, 101, 101, 101]
styles = ['r--','-','-','-','-']
for file, label, point, style in zip(files,labels, norm_points, styles):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,-1] = smooth_data(data[:,-1], 100)
   max_exper = get_scale_factor(exper, energy_value=point)
   scale_factor = max_exper/get_scale_factor(data, energy_value=point, col=-1)
   plt.semilogy(data[:,0], scale_factor*data[:,-1], style, linewidth='1.5', label=label)

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')
plt.tight_layout()


plt.figure(4)
# plt.title('Figure 3 from paper')
plt.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data")

files=[
         'Imin_95/full_area_0.05_step_1.csv',
         'Imin_95/full_area_0.05_step_2.csv',
         'Imin_95/full_area_0.05_step_5.csv',
         'Imin_95/full_area_0.05_step_10.csv',
         'Imin_95/full_area_0.05_step_12.csv'
         ]


labels=[
         '$\Delta$I = 1x10$^11$ Wcm$^{-2}$',
         '$\Delta$I = 2x10$^11$ Wcm$^{-2}$',
         '$\Delta$I = 5x10$^11$ Wcm$^{-2}$',
         '$\Delta$I = 10x10$^11$ Wcm$^{-2}$',
         '$\Delta$I = 12x10$^11$ Wcm$^{-2}$',
         ]

norm_points=[101, 101, 101, 101, 101]
styles = ['r--','-','-','-','-']
for file, label, point, style in zip(files,labels, norm_points, styles):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,-1] = smooth_data(data[:,-1], 100)
   max_exper = get_scale_factor(exper, energy_value=point)
   scale_factor = max_exper/get_scale_factor(data, energy_value=point, col=-1)
   plt.semilogy(data[:,0], scale_factor*data[:,-1], style, linewidth='1.5', label=label)

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')
plt.tight_layout()

# ============================================================

plt.figure(5)
# plt.title('Plot Figure 4 from paper')
plt.semilogy(exper[:,0],exper[:,1], 'k', linewidth='1.5', label="Experimental Data")

files=[
         'Imin_95/full_area_0.01_step_1.csv',
         'Imin_95/full_area_0.03_step_1.csv',
         'Imin_95/full_area_0.05_step_1.csv',
         'Imin_95/full_area_0.07_step_1.csv',
         'Imin_95/full_area_0.09_step_1.csv',
         'Imin_95/full_area_0.1_step_1.csv'

         ]


labels=['0.01mm$^2$',
         '0.03mm$^2$',
         '0.05mm$^2$',
         '0.07mm$^2$',
         '0.09mm$^2$', 
         '0.1mm$^2$' 
         ]

norm_points= np.ones(len(labels))*101
for file, label, point in zip(files,labels, norm_points):
   data = np.loadtxt(file, skiprows=1, delimiter=',')
   data[:,-2] = smooth_data(data[:,-2])
   max_exper = get_scale_factor(exper, energy_value=point)
   scale_factor = max_exper/get_scale_factor(data, energy_value=point, col=-2)
   plt.semilogy(data[:,0], scale_factor*data[:,-2], linewidth='1.5', label=label)

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-2])
plt.legend(fontsize='small')
plt.tight_layout()


plt.show()