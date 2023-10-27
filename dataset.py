import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import pandas as pd

#extracting useful features
def og_evaluate(filename):
    hdulist = fits.open(filename)
    times, rates = hdulist[1].data["TIME"], hdulist[1].data["RATE"]
    hdulist.close()

    found = []
    g = gaussian_filter(rates, sigma=10)
    peaks, _ = find_peaks(g)
                
    prominences, _, _ = peak_prominences(g, peaks)
    if len(prominences) == 0:
        return []
                
    selected = prominences > 0.5 * (np.min(prominences) + np.max(prominences))
    if len(selected) == 0:
        return []
    
    top = peaks[selected]
    if len(top) > 10:
        return []
    eigth_peak_widths = peak_widths(g, top, rel_height=0.8)        
                
    per = np.percentile(g, 75)

    for i in range(len(g[top])):
        found.append({
            'start_time': times[int(eigth_peak_widths[2][i])],
            'end_time': times[int(eigth_peak_widths[3][i])],
            'max_time': times[top[i]],
            'peak_width': eigth_peak_widths[0][i],
            'peak_risetime': top[i] - eigth_peak_widths[2][i],
            'peak_falltime': eigth_peak_widths[3][i] - top[i],
            'peak_flux': rates[top[i]],
        })
    
    return found

#First we downloaded all the .lc files and created a text file containing the paths. Now this code runs through every file and creates a dataframe with above features
def process_all_files(file_paths_file):
    with open(file_paths_file, 'r') as file:
        file_paths = file.read().splitlines()

    df = pd.DataFrame()

    for file_path in file_paths:
        found = og_evaluate(file_path)
        df = df.append(found, ignore_index=True)

    df.to_csv('feature_data.csv', index=False)

file_paths_file = 'lc_file_paths.txt'

process_all_files(file_paths_file)
