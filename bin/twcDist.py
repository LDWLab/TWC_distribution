#!/usr/bin/env python3
import sys, csv
import numpy as np
from twincons.TwinCons import data_to_diverging_gradients
from scipy.cluster.vq import kmeans
import scipy.cluster as cluster
import matplotlib.pyplot as plt

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

def calculate_stats(twcListData, ks):
    centroids, avg_distance = kmeans(twcListData, ks)
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
    return [np.max(lowOutlierDatas), groups]

def plot_data_dist_and_anomaly_threshold(data, lower_limit, lowerStd, outpath):
    plt.title(f"Threshold at {round(lower_limit,2)} determined with 7 kmeans.")
    plt.hist(list(data.values()), bins=50)
    plt.axvline(x=lower_limit, color='r')
    plt.axvspan(lower_limit-lowerStd, lower_limit+lowerStd, color='r', alpha=.5)
    plt.savefig(outpath)
    plt.close()
    return True

def main(inTWC, outpng, outcsv):
    twc_data = read_twc_data(inTWC)
    lower_limit, upper_limit, twc_anomally = find_anomalies(twc_data, 1)
    
    twcListData = list(twc_data.values())

    i = 0
    statMins = list()
    while i < 100:
        statMins.append(calculate_stats(twcListData, 7)[0])
        i += 1

    groups = calculate_stats(twcListData, 7)[1]
    print(round(np.mean(statMins),3), round(np.std(statMins),3))
    signatureThreshold = np.mean(statMins)

    plt.title(f"TWC signature threshold {round(signatureThreshold,2)}")
    plt.scatter(twcListData, np.arange(0,len(twcListData)), c=groups)
    plt.axvline(x=signatureThreshold, c='r')
    plt.axvspan(signatureThreshold-np.std(statMins), signatureThreshold+np.std(statMins), color='r', alpha=.5)
    plt.savefig(outpng.replace('.png', '_groups.png'))
    plt.close()

    plot_data_dist_and_anomaly_threshold(twc_data, signatureThreshold, np.std(statMins), outpng)
    with open(outcsv, 'w') as f:
        f.write(f"TWC signature threshold and STD {round(np.mean(statMins),3)},{round(np.std(statMins),3)}\n")
        f.write("Residue,TWC\n")
        for key in twc_anomally.keys():
            f.write(f"{key},{twc_anomally[key]}\n")


main(sys.argv[1], sys.argv[2], sys.argv[3])