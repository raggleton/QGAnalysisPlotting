#!/usr/bin/env python


"""
Figure out binning factors to get N bins between certain bin edges,
such that the bin width is guaranteed to be a certain fraction of the bin edge
"""

from math import log

edges = [
50,
65,
88,
109,
180,
250,
318,
388,
479,
535,
]

# fractional bin width cannot be below this
# lambda so you can do it as a function of pt
min_multiplier = lambda x: 1.25

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
    print("multiplier: %.3f" % smallest_multiplier, "gives", N_bins, "bins:", this_bin_edges, compromise_str)
