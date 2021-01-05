#!/usr/bin/env python3
import sys, csv
import numpy as np
from twincons.TwinCons import data_to_diverging_gradients
from scipy.cluster.vq import kmeans, kmeans2
import scipy.cluster as cluster
import matplotlib.pyplot as plt

def plot_data_dist_and_anomaly_threshold(data, lower_limit, upper_limit, std, outpath):
    plt.title(f"Thresholds determined with {std} standard deviations from the mean.")
    plt.hist(list(data.values()), bins=50)
    plt.axvline(x=lower_limit, color='r')
    plt.axvline(x=upper_limit, color='g')
    plt.savefig(outpath)
    plt.close()
    return True

def find_anomalies(data, std_devs = 2):
    
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

    return lower_limit, upper_limit, anomalies

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

def read_twc_data(twc_loc):
    twc_data = dict()
    with open(twc_loc) as twc_f:
        reader = csv.reader(twc_f)
        for line in reader:
            if line[0] == 'resNum':
                continue
            if line[1] == 'NA':
                continue
            twc_data[line[0]] = float(line[1])
    return twc_data

def main(inTWC, outpng, outcsv):
    twc_data = read_twc_data(inTWC)
    lower_limit, upper_limit, twc_anomally = find_anomalies(twc_data, 1)
    
    twcListData = list(twc_data.values())
    centroids, avg_distance = kmeans(twcListData, 7)
    groups, cdist = cluster.vq.vq(twcListData, centroids)

    minGroup = groups[np.argmin(twcListData)]
    maxGroup = groups[np.argmax(twcListData)]

    minLoc, maxLoc = list(), list()
    for ix, g in enumerate(groups):
        if g == minGroup:
            minLoc.append(ix)
        if g == maxGroup:
            maxLoc.append(ix)

    lowOutlierDatas = [twcListData[loc] for loc in minLoc]
    return np.max(lowOutlierDatas)
    #plt.scatter(twcListData, np.arange(0,len(twcListData)), c=groups)
    #plt.savefig('test2.png')
    #plt.close()

    #plot_data_dist_and_anomaly_threshold(twc_data, lower_limit, upper_limit, 1, outpng)
    #with open(outcsv, 'w') as f:
    #    f.write("Residue,TWC\n")
    #    for key in twc_anomally.keys():
    #        f.write(f"{key},{twc_anomally[key]}\n")

i = 0
statMins = list()
while i < 100:
    statMins.append(main(sys.argv[1], sys.argv[2], sys.argv[3]))
    i += 1

print(np.mean(statMins), np.std(statMins))


