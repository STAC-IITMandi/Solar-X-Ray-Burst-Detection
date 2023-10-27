import astropy
import scipy
from astropy.io import fits
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from astropy.table import Table
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_prominences

class x_ray_burst:
    def __init__(self, gti_filename, lc_filename):
        self.gti_filename = gti_filename
        self.lc_filename = lc_filename
    def data(self):
        gti_hdulist = fits.open(gti_filename)
        gti_data = gti_hdulist[1].data
        gti_start_times = gti_data['START']  
        gti_end_times = gti_data['STOP']  
    
        lc_hdulist = fits.open(lc_filename)
        lc_data = lc_hdulist[1].data
        lc_times = lc_data['TIME']  
        lc_rates = lc_data['RATE']  
        return gti_start_times, gti_end_times, lc_times, lc_rates
    def combined_data(self):
        gti_start_times, gti_end_times, lc_times, lc_rates = self.data()
        filtered_lc_times = []
        filtered_lc_values = []
        for i in range(len(lc_times)):
            for j in range(len(gti_start_times)):
                if gti_start_times[j] <= lc_times[i] <= gti_end_times[j]:
                    filtered_lc_times.append(lc_times[i])
                    filtered_lc_values.append(lc_rates[i])
                    break
        return filtered_lc_times, filtered_lc_values
    def plot_combined_data(self):
        times, values = self.combined_data()
        plt.plot(times, values)
        plt.xlabel("Time")
        plt.ylabel("XRay Flux")
        plt.show()
    def peaks(self):
        times, values = self.combined_data()
        smoothed_values = gaussian_filter1d(values, sigma=5)  
        plt.figure(figsize=(12, 6))
        plt.plot(times, smoothed_values, label='Smoothed Data (Gaussian)')
        plt.xlabel("Time")
        plt.ylabel("XRay Flux")
        plt.legend()
        plt.grid(True)
        plt.show()
        peaks, _= find_peaks(smoothed_values, height=10, prominence=1000)
        prominences= peak_prominences(smoothed_values, peaks)[0]
        plt.figure(figsize=(10,6))
        plt.plot(smoothed_values, label= 'data')
        plt.plot(peaks, smoothed_values[peaks], 'ro', label= 'detected peaks')
        plt.xlabel('time/ data point index')
        plt.ylabel('value')
        plt.legend()
        plt.grid(True)
        plt.show()
    def features(self):
        times, smoothed_values= self.combined_data()
        timestamps = np.array(times)
        flux = np.array(smoothed_values)
        threshold = 1
        burst_start = np.where(flux > threshold)[0][0]
        burst_end = np.where(flux > threshold)[0][-1]
        burst_duration = timestamps[burst_end] - timestamps[burst_start]

        rise_threshold = 0.5
        decay_threshold = 0.5

        peak_flux = np.max(flux)
        time_of_rise= timestamps[np.where(flux>rise_threshold)[0][0]]
        time_of_decay= timestamps[np.where(flux>decay_threshold)[0][-1]]

        print(f"duration of burst: {burst_duration}")
        print(f"time of rise: {time_of_rise}")
        print(f'time of decay: {time_of_decay}')
        print(f"peak flux: {peak_flux}")

    
x= input("Enter the GTI file path: ")
y= input("Enter the LC file path: ")
gti_filename= x.strip('"')
lc_filename= y.strip('"')

burst = x_ray_burst(gti_filename,lc_filename)
burst.plot_combined_data()
burst.peaks()
times, smoothed_values = burst.combined_data()
burst.features()