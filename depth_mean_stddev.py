#!/usr/bin/env python3

"""For a mosdepth bed file provided to stdin, determine the mean and stddev as
if the data were normal, and report these values.
Usage:
    tabix ${BAM%.*}.median.regions.bed.gz <list autosomes> | depth_mean_stddev.py
Example:
    tabix ${BAM%.*}.median.regions.bed.gz {1..22} | depth_mean_stddev.py
"""

import sys

import numpy as np
from scipy.stats import mode, norm

# pull coverage depth into an array
depth = np.loadtxt(sys.stdin, delimiter="\t", usecols=3, dtype=np.float64)

# calculate mode of values > 0
nonzero_mode = mode(depth[depth > 0])[0]

# calculate mu, sigma of values from 1 to 2*mode
mean, stddev = norm.fit(depth[(depth > 0) & (depth < (2 * nonzero_mode))])

# output mu, sigma
print(mean, stddev)
