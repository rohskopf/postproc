"""
Calculate power spectrum of mode coordinates.

Loading with np.loadtxt takes 11 minutes to load xm.dat.
Loading with pd takes 4 minutes to load xm.dat.
"""
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


"""
Declare settings.
"""
#dirname = "e_mode_100_1/xm.dat"
dirname = "100K"
# At what intervals time points are sampled
samplingInterval = 2.5e-3 # data collected every 25 fs = 2.5e-3 ps
# How many time points there are, i.e. sampling frequency
samplingFrequency = 1/samplingInterval

Path("ps_"+dirname).mkdir(parents=True, exist_ok=True)

"""
Read the file.
"""
start_time = time.time()
#xm = np.loadtxt("e_mode_100_1/xm.dat")
data_xm = pd.read_csv(dirname+"/xm.dat", header=None, delim_whitespace=True).to_numpy()
#data_vm = pd.read_csv(dirname+"/vm.dat", header=None, delim_whitespace=True).to_numpy()

execution_time = (time.time() - start_time)
execution_time = execution_time/60. # minutes
print("Took %f minutes to load xm." % (execution_time))

# Begin time period of the signals
beginTime = data_xm[0,0];
# End time period of the signals
data_shape = np.shape(data_xm)
npoints = data_shape[0] #-350000
nmodes = data_shape[1]-1
print("Number of modes (columns): %d" % (nmodes))
endTime = data_xm[npoints-1,0]
print("end time: %f" % (endTime))
# Time points
time_points = np.arange(beginTime, endTime, samplingInterval);

"""
Load the indices of interest.
"""
indices = np.loadtxt("INDICES_MODECOOR_PS").astype(int)
#print(indices)

"""
Loop through indices and calculate PS 
"""
start_time = time.time()
for indx in indices:

  print(indx)

  # Create subplot
  figure, axis = plt.subplots(1, 1)
  plt.subplots_adjust(hspace=1)

  """
  Mode amplitude
  """
  # Get the signal
  amplitude = data_xm[0:npoints-1,indx+1] # Need to index +1 beacuse first column is time.
  mean_amp = np.mean(amplitude)
  amplitude = amplitude - mean_amp
  # Frequency domain representation
  fourierTransform = np.fft.fft(amplitude)/len(amplitude)           # Normalize amplitude
  fourierTransform = fourierTransform[range(int(len(amplitude)/2))] # Exclude sampling frequency
  tpCount     = len(amplitude)
  values      = np.arange(int(tpCount/2))
  timePeriod  = tpCount/samplingFrequency
  frequencies = values/timePeriod
  # Axes
  #axis[0].plot(time_points, amplitude+mean_amp)
  #axis[0].set_xlabel('Time')
  #axis[0].set_ylabel('X sqrt(M)*A')
  axis.plot(frequencies, abs(fourierTransform)**2)
  axis.set_xlabel('Frequency')
  axis.set_ylabel('Amplitude')
  axis.set_xlim(0,30)
  """
  Mode velocity.
  """
  
  """
  # Get the signal
  amplitude = data_vm[0:npoints-1,indx+1] # Need to index +1 beacuse first column is time.
  mean_amp = np.mean(amplitude)
  amplitude = amplitude - mean_amp
  # Frequency domain representation
  fourierTransform = np.fft.fft(amplitude)/len(amplitude)           # Normalize amplitude
  fourierTransform = fourierTransform[range(int(len(amplitude)/2))] # Exclude sampling frequency
  #tpCount     = len(amplitude)
  #values      = np.arange(int(tpCount/2))
  #timePeriod  = tpCount/samplingFrequency
  #frequencies = values/timePeriod
  # Create subplot
  #figure, axis = plt.subplots(2, 1)
  #plt.subplots_adjust(hspace=1)
  # Axes
  axis[2].plot(time_points, amplitude+mean_amp)
  axis[2].set_xlabel('Time')
  axis[2].set_ylabel('V sqrt(M)*A/ps')
  axis[3].plot(frequencies, abs(fourierTransform)**2)
  axis[3].set_xlabel('Frequency')
  axis[3].set_ylabel('Amplitude')
  axis[3].set_xlim(0,30)
  """
  
  figname = "ps_"+dirname+"/mode%d.png" % (indx)
  plt.savefig(figname)

execution_time = (time.time() - start_time)
execution_time = execution_time/60. # minutes
print("Took %f minutes to calculate power spectrums." % (execution_time))
