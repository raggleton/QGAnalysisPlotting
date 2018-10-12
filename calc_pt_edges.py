#!/usr/bin/env python


"""
Figure out binning factors to get N bins between certain bin edges,
such that the bin width is guaranteed to be a certain fraction of the bin edge
"""

from math import log

# AK4
print("Doing AK4")
edges = [
88.0,
114.0,
185.0,
254.0,
320.0,
390.0,
481.0,
539.0,
]

# AK8
# print("Doing AK8")
# edges = [
# 120.0,
# 148.0,
# 233.0,
# 326.0,
# 405.0,
# 477.0,
# # 574.0,
# 614.0,
# ]

# print("Doing AK8 specialised")
# edges = [
# 68,
# 88,
# 129,
# 185,
# 252,
# 320,
# 390,
# 479,
# 536
# ]

# fractional bin width cannot be below this
# lambda so you can do it as a function of pt
min_multiplier = lambda x: 1.25

all_bins = []

for low, high in zip(edges[:-1], edges[1:]):
    smallest_multiplier = 9999
    N_bins = 0
    for N in range(10, 0, -1):
        a = pow(high/low, 1./N)
        if a >= min_multiplier(low) and a < smallest_multiplier:
            smallest_multiplier = a
            N_bins = N

    compromise_str = ""
    if smallest_multiplier > 100:
        smallest_multiplier = min_multiplier(low)
        N_bins = 1
        compromise_str = "** COMPROMISE, requires multiplier = %.3f" % (high/low)

    this_bin_edges = [round(low*pow(smallest_multiplier, i)) for i in range(0, N_bins+1)]
    all_bins += this_bin_edges
    print("multiplier: %.3f" % smallest_multiplier, "gives", N_bins, "bins:", this_bin_edges, compromise_str)

while all_bins[-1] < 10000:
    all_bins.append(round(all_bins[-1] * 1.3))

# set cos I'm lazy and didnt filter out repeats
print("All bins:", sorted(list(set(all_bins))))