#!/usr/bin/env python3
import sys, csv
import numpy as np
#from twincons.TwinCons import data_to_diverging_gradients
from scipy.cluster.vq import kmeans
import scipy.cluster as cluster
import matplotlib
import matplotlib.pyplot as plt

def data_to_diverging_gradients(datapoint, maxdata, mindata, positivegradient, negativegradient):
    '''Maps a data point to a diverging colormap depending on whether its above or bellow 0.
    Returns a hex code.
    '''
    if datapoint == 'NA':
        return '#808080'
    if datapoint >= 0:
        grad = np.atleast_2d(np.linspace(0,datapoint/maxdata,256)).T
        rgb = plt.get_cmap(positivegradient)(grad[len(grad)-1])[0][:3]
    else:
        grad = np.atleast_2d(np.linspace(0,datapoint/mindata,256)).T
        rgb = plt.get_cmap(negativegradient)(grad[len(grad)-1])[0][:3]
    return matplotlib.colors.rgb2hex(rgb)

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
    highOutlierDatas = [twcListData[loc] for loc in maxLoc]
    return [np.max(lowOutlierDatas), np.min(highOutlierDatas), groups]

def plot_data_dist_and_anomaly_threshold(data, lowerLimit, lowerStd, upperlimit, upperStd, outpath):
    plt.title(f"Signature threshold at {round(lowerLimit,2)} determined with 5 kmeans.\nConserved threshold at {round(upperlimit,2)}")
    plt.hist(list(data.values()), bins=50)
    plt.axvline(x=lowerLimit, color='r')
    plt.axvline(x=upperlimit, color='g')
    plt.axvspan(lowerLimit-lowerStd, lowerLimit+lowerStd, color='r', alpha=.3)
    plt.axvspan(upperlimit-upperStd, upperlimit+upperStd, color='g', alpha=.3)
    plt.savefig(outpath)
    plt.close()
    return True

def main(inTWC, outpng, outcsv, numkmeans):
    twc_data = read_twc_data(inTWC)
    lower_limit, upper_limit, twc_anomally = find_anomalies(twc_data, 1)
    numkmeans = int(numkmeans)
    twcListData = list(twc_data.values())

    i = 0
    statMins, statsMaxs = list(), list()
    while i < 100:
        statMins.append(calculate_stats(twcListData, numkmeans)[0])
        statsMaxs.append(calculate_stats(twcListData, numkmeans)[1])
        i += 1

    groups = calculate_stats(twcListData, numkmeans)[2]
    print(f"Min thr: {round(np.mean(statMins),3)}, STD: {round(np.std(statMins),3)}\tMax thr: {round(np.mean(statsMaxs),3)}, STD: {round(np.std(statsMaxs),3)}")
    signatureThreshold = np.mean(statMins)
    conservedThreshold = np.mean(statsMaxs)

    plt.title(f"TWC signature threshold {round(signatureThreshold,2)}")
    plt.scatter(twcListData, np.arange(0,len(twcListData)), c=groups)
    plt.axvline(x=signatureThreshold, c='r')
    plt.axvline(x=conservedThreshold, c='g')
    plt.axvspan(signatureThreshold-np.std(statMins), signatureThreshold+np.std(statMins), color='r', alpha=.3)
    plt.axvspan(conservedThreshold-np.std(statsMaxs), conservedThreshold+np.std(statsMaxs), color='g', alpha=.3)
    plt.savefig(outpng.replace('.png', '_groups.png'))
    plt.close()

    plot_data_dist_and_anomaly_threshold(twc_data, signatureThreshold, np.std(statMins), conservedThreshold, np.std(statsMaxs), outpng)
    with open(outcsv, 'w') as f:
        f.write(f"TWC signature threshold and STD {round(np.mean(statMins),3)},{round(np.std(statMins),3)}\n")
        f.write("Residue,TWC\n")
        for key in twc_anomally.keys():
            f.write(f"{key},{twc_anomally[key]}\n")


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])