#!/usr/bin/env python3

# (C) 2015, Alexander Korsunsky

import numpy as np

import scipy.spatial.distance
import scipy.cluster.hierarchy
import scipy.cluster.vq

import matplotlib.pyplot as plt


def fpkm_filter(array, threshold):
    outarray = array.copy()
    outarray[outarray < threshold] = 0  # zero out all values below a threshold
    outarray = outarray[~np.all(outarray == 0, axis=1)]  #  weed out zero lines
    return outarray


def calculate_linkage(observations, distance_method, linkage_method):
    linkage = scipy.cluster.hierarchy.linkage(
        scipy.spatial.distance.pdist(observations, distance_method), linkage_method)
    return linkage


def load_fpkm(cuffcmp_tracking_file, labelfile, use_labels=None):
    labels = np.genfromtxt(labelfile, dtype=str)

    if use_labels:
        used_indices = np.in1d(labels, use_labels)
    else:
        used_indices = np.ones(len(labels), dtype=bool)

    fpkmmatrix = np.load(cuffcmp_tracking_file)


    return fpkmmatrix[used_indices], labels[used_indices]

def plot_hierarchical(cuffcmp_tracking_file, labelfile, distance_method, linkage_method, use_labels=None):
    fpkmmatrix, labels = load_fpkm(cuffcmp_tracking_file, labelfile, use_labels)
    fpkmmatrix = fpkm_filter(fpkmmatrix, 50.0)

    Z = calculate_linkage(fpkmmatrix.T, distance_method, linkage_method)

    dend = scipy.cluster.hierarchy.dendrogram(Z, labels=labels, orientation='right', leaf_rotation=45,
                                              distance_sort=True, truncate_mode='lastp', p=26)


    # matplotlib.pyplot.figure(6)

def do_kmeans(cuffcmp_tracking_file, labelfile, num_clusters, use_labels=None):
    fpkmmatrix, labels = load_fpkm(cuffcmp_tracking_file, labelfile, use_labels)
    fpkmmatrix = fpkm_filter(fpkmmatrix, 50.0)


    centroids,_ = scipy.cluster.vq.kmeans2(fpkmmatrix.T, num_clusters, minit='points')
    idx,_ = scipy.cluster.vq.vq(fpkmmatrix.T,centroids)

    for k in range(0,num_clusters):
        print("\n", k);
        print('\n'.join(labels[idx == k]))


# idx_0m = np.array([9,18,21,22,28,31,32,37,38,41,47,54,59])-1
# idx_4m = np.array([7,8,13,17,25,27,30,33,36,40,42,49,50])-1
# idx_5m = np.array([6])-1
# idx_4m = np.array([7,8,13,17,27,30,33,36,40,42,50])-1
# idx_8m = np.array([5,11,12,16,19,20,29,43,44,48])-1
# idx_12m = np.array([3,4,10,14,26,34,35,39,45,51,52,56])-1
# idx_12m = np.setdiff1d(idx_12m,idx_controls)
# idx_14m = np.array([25,49])-1
# idx_15m = np.array([1,46])-1
#
# idx_negatives = np.union1d(idx_controls,idx_0m)
# idx_4m5m = np.union1d(idx_4m,idx_5m)
# idx_14m5m = np.union1d(idx_4m,idx_5m)
#
#
# fpkmmatrix = np.load("/export/bse/diff/cuffcompare-all-bg/cuffcmp.tracking.fpkmmatrix.npy")
# binmatrix = np.array(fpkmmatrix, dtype=bool)
#
# bin_negatives = binmatrix[:,idx_negatives]
# fpkmmatrix_negatives_removed = fpkmmatrix[~np.any(bin_negatives, axis=1),:
#
#
# fpkm_4m5m = fpkmmatrix_negatives_removed[:,idx_4m5m]
# fpkm_8m = fpkmmatrix_negatives_removed[:,idx_8m]
# fpkm_12m = fpkmmatrix_negatives_removed[:,idx_12m]
# fpkm_14m15m = fpkmmatrix_negatives_removed[:,idx_14m15m]
