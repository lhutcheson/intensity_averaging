import numpy as np
import matplotlib.pyplot as plt
import glob

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

def strip_num(filename):
    import re
    pat=re.compile("[0-9]{1}.[0-9]{3}")
    m=pat.search(filename)
    return(float(m.group(0)))

naive_list = glob.glob("full_coherent_area_0.1_step_*")
naive_list.sort()
full_list = glob.glob("naive_coherent_area_0.1_step_*")
full_list.sort()

exper = np.loadtxt('experimental_data.csv', skiprows=1, delimiter=',')
max_exper = get_scale_factor(exper, 90)
plt.semilogy(exper[:,0], exper[:,1], 'k', label="Experimental Data")

colours=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

for count, file in enumerate(naive_list):
   step=strip_num(file)
   # data = np.loadtxt(file, skiprows=1, delimiter=',')
   # if args.smooth_data:
   smooth = np.loadtxt(file, skiprows=1, delimiter=',')
   smooth[:,1] = smooth_data(smooth[:,1])
   # max_data = get_scale_factor(data, 93.2)
   max_smooth = get_scale_factor(smooth, 104.5)

   # scale_factor = max_exper/max_data
   # plt.semilogy(data[:,0], scale_factor*data[:,1], alpha=0.5, color=colours[count])

   scale_factor = max_exper/max_smooth
   plt.semilogy(smooth[:,0], scale_factor*smooth[:,1], label='$\Delta I = $'+str(step)+'$x10^{14}$', color=colours[count])
   # else:
   #    max_data = get_scale_factor(data, 90)
   #    scale_factor = max_exper/max_data
   #    plt.semilogy(data[:,0], scale_factor*data[:,1], label='$\Delta I $= '+str(step)+'$x10^{14}$')

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-1])
plt.legend(fontsize='small')

plt.show()

plt.semilogy(exper[:,0], exper[:,1], 'k', label="Experimental Data")

scale = [101.3, 100, 98.2, 96.6, 93.2, 93.3]
for count, file in enumerate(full_list):
   step=strip_num(file)
   # data = np.loadtxt(file, skiprows=1, delimiter=',')
   # if args.smooth_data:
   smooth = np.loadtxt(file, skiprows=1, delimiter=',')
   smooth[:,1] = smooth_data(smooth[:,1])
   # max_data = get_scale_factor(data, scale[count])
   max_smooth = get_scale_factor(smooth, scale[count])

   # scale_factor = max_exper/max_data
   # plt.semilogy(data[:,0], scale_factor*data[:,1], alpha=0.5, color=colours[count])

   scale_factor = max_exper/max_smooth
   plt.semilogy(smooth[:,0], scale_factor*smooth[:,1], label='$\Delta I = $'+str(step)+'$x10^{14}$', color=colours[count])
   # else:
   #    max_data = get_scale_factor(data, 90)
   #    scale_factor = max_exper/max_data
   #    plt.semilogy(data[:,0], scale_factor*data[:,1], label='$\Delta I $= '+str(step)+'$x10^{14}$')

plt.xlabel('Photon Energy (eV)')
plt.ylabel('Harmonic Intensity (arb. units)')
plt.xlim([0,160])
plt.ylim([10e-7, 10e-1])
plt.legend(fontsize='small')

plt.show()