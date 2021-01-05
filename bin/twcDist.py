#!/usr/bin/env python3
import sys, csv
import numpy as np
from twincons.TwinCons import data_to_diverging_gradients


def plot_data_dist_and_anomaly_threshold(data, lower_limit, upper_limit, std, outpath):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    plt.title(f"Thresholds determined with {std} standard deviations from the mean.")
    plt.hist(list(data.values()), bins=50)
    plt.axvline(x=lower_limit, color='r')
    plt.axvline(x=upper_limit, color='g')
    plt.savefig(outpath)
    plt.close()
    return True

def find_anomalies(data, outpath, std_devs = 2):
    
    anomalies = dict()
    data_list = list(data.values())
    random_data_std = np.std(data_list)
    random_data_mean = np.mean(data_list)
    anomaly_cut_off = random_data_std * std_devs
    lower_limit  = random_data_mean - anomaly_cut_off
    upper_limit = random_data_mean + anomaly_cut_off

    for pdbRef, outlier in data.items():
        if outlier > upper_limit or outlier < lower_limit:
            anomalies[pdbRef] = outlier
    
    plot_data_dist_and_anomaly_threshold(data, lower_limit, upper_limit, std_devs, outpath)

    return anomalies

def remove_anomalies(subfspec, anomalies):
    clean_subfspec = dict()
    for key, value in subfspec.items():
        if key not in anomalies.keys():
            clean_subfspec[key] = value
    return clean_subfspec

def add_highly_conserved_and_anomalies(subf_clean, minval, addition_dict):
    for k in addition_dict.keys():
        if k not in subf_clean.keys():
            subf_clean[k] = minval
    return subf_clean

def reflect_dict_data(data_dict):
    reflect_dict = dict()
    for k, v in data_dict.items():
        reflect_dict[k] = v*-1
    return reflect_dict

def data_to_colors(final_data_dict):
    ze_max, ze_min = max(final_data_dict.values()), min(final_data_dict.values())
    color_dict = dict()
    for k, v in final_data_dict.items():
        color_dict[k] = data_to_diverging_gradients(v, ze_max, ze_min, 'Greens', 'Purples')
    return color_dict

def twc_fun(twc_loc, outpath):
    import csv
    twc_data = dict()
    with open(twc_loc) as twc_f:
        reader = csv.reader(twc_f)
        for line in reader:
            if line[0] == 'resNum':
                continue
            if line[1] == 'NA':
                continue
            twc_data[line[0]] = float(line[1])
    twc_anomally = find_anomalies(twc_data, outpath, 1)
    return twc_anomally


twc_anomally = twc_fun(sys.argv[1], sys.argv[2])
with open(sys.argv[3], 'w') as f:
    f.write("Residue,TWC\n")
    for key in twc_anomally.keys():
        f.write(f"{key},{twc_anomally[key]}\n")
