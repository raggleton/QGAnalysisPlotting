#!/usr/bin/env python

"""Plot the jets/pf and genjets/genparticles for good & bad herwig events"""

import matplotlib.pyplot as plt
import numpy as np

from pf import output_bad as pf_output
from gp import output_bad as gp_output
from genjet import output_bad as genjet_output
from jet import output_bad as jet_output

# from pf import output_good as pf_output
# from gp import output_good as gp_output
# from genjet import output_good as genjet_output
# from jet import output_good as jet_output


def create_particle_array(raw_data):
    particles = []
    for line in raw_data.splitlines():
        parts = [float(x.strip()) for x in line.split(":")[1:4]]
        particles.append(parts)

    particles = np.array(particles)

    if particles.shape[1] != 3:
        raise IndexError("Each entry doesn't have 3 components")

    return particles


pf_particles = create_particle_array(pf_output)
gp_particles = create_particle_array(gp_output)
genjets = create_particle_array(genjet_output)
jets = create_particle_array(jet_output)

# filtering

jet_pt_cut = 30
jets = jets[jets[:, 0] > jet_pt_cut]

print pf_particles.shape[0], ' PF particles'
print gp_particles.shape[0], ' GenParticles'
print genjets.shape[0], ' GenJets'
print jets.shape[0], ' Jets'

# GEN
plt.scatter(gp_particles[:,1], gp_particles[:,2], marker='o', label='GenParticles', alpha=0.4, s=8*gp_particles[:,0], c='r')
plt.scatter(genjets[:,1], genjets[:,2], marker='x', label='GenJets', c='r', s=2*genjets[:,0])
# RECO
plt.scatter(pf_particles[:,1], pf_particles[:,2], marker='s', label='PF particles', alpha=0.4, s=2*pf_particles[:,0], c='b')
plt.scatter(jets[:,1], jets[:,2], marker='x', label='Jets (pT > %d)' % jet_pt_cut, c='b', s=2*jets[:,0])

plt.gca().tick_params(which='both', bottom=True, top=True)
plt.grid(which='both', linestyle='dotted')
plt.xlabel('eta')
plt.ylabel('phi')
plt.xlim(-3, 3)
plt.ylim(-np.pi, np.pi)
plt.legend()

plt.suptitle('Herwig QCD, bad, event: 9719071, lumi: 67493')
plt.savefig("plot_bad.pdf")

# plt.suptitle('Herwig QCD, good, event: 52710, lumi: 367')
# plt.savefig("plot_good.pdf")
