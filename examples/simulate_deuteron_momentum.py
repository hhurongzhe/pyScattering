# a deuteron solver in momentum space.
# for soft interactions (em and emn in momentum space), kmax~1000MeV would be enough;
# for hard interactions (idaho and av18 in position space, cdbonn), kmax~4000MeV would be enough.


import sys

sys.path.append("./lib")
import utility
import profiler
import time

import numpy as np
import nn_studio as nn_studio
import chiral_potential as chiral_potential
from scipy import linalg
from scipy.special import spherical_jn as jn
import matplotlib.pyplot as plt
import constants as const

utility.header_message()

################################################################################################################

utility.section_message("Initialization")

# initialize an object for computing T-matrices, phase shifts,
nn = nn_studio.nn_studio(jmin=0, jmax=1, tz=0, Np=130, mesh_type="gauleg_finite")

potential_type = "n4lo+sms450"

print(f"potential type : {potential_type}\n")
# initialize an object for the chiral interaction
potential = chiral_potential.two_nucleon_potential(potential_type)

# give the potential to the nn-analyzer
nn.V = potential

################################################################################################################


mu = const.Mp * const.Mn / (const.Mp + const.Mn)

N = 2 * (nn.Np)

H = np.zeros((N, N))
T = np.zeros((N, N))

ww = np.hstack((nn.wmesh, nn.wmesh))
pp = np.hstack((nn.pmesh, nn.pmesh))

print("Paramaters for solving deuteron\n")
print(f"k_min (MeV) : {nn.pmesh[0]}")
print(f"k_max (MeV) : {nn.pmesh[-1]}")
print(f"k mesh number : {nn.Np}")

################################################################################################################


################################################################################################################
utility.section_message("Building Potential")

t1 = time.time()
is_coupled, channel = True, (1, 0, 0)
V = nn.setup_Vmtx(is_coupled, channel)
t2 = time.time()
profiler.add_timing("Set up V", t2 - t1)
################################################################################################################


################################################################################################################
utility.section_message("Building Hamiltonian")

for i, p_bra in enumerate(pp):
    for j, p_ket in enumerate(pp):
        Tij = 0
        if i == j:
            Tij = p_bra**2 / (2 * mu)
            T[i][j] = Tij

        V[i][j] = V[i][j] * p_bra * p_ket * np.sqrt(ww[i] * ww[j])
H = T + V
t3 = time.time()
profiler.add_timing("Set up H", t3 - t2)
################################################################################################################


################################################################################################################
utility.section_message("Solving Schrodinger")

eigvals, eigvecs = linalg.eigh(H)
t4 = time.time()
profiler.add_timing("Diagonalization", t4 - t3)
################################################################################################################


################################################################################################################
utility.section_message("Calculating Observables")

s = np.argsort(eigvals)
E = eigvals[s[0]]
psi_k = eigvecs[:, s[0]]

print(f"E = {E}")


def phi_r(this_r):
    psi_r_s = 0
    psi_r_d = 0
    for i, p in enumerate(pp):
        if i < N / 2:
            j = jn(0, this_r * p / const.hbarc)
            psi_r_s += np.sqrt(ww[i]) * p * j * psi_k[i]
        else:
            j = jn(2, this_r * p / const.hbarc)
            psi_r_d += np.sqrt(ww[i]) * p * j * psi_k[i]

    return psi_r_s, psi_r_d


rmin = 1e-16
rmax = 40
rr, wr = nn.gauss_legendre_line_mesh(rmin, rmax)
psi_r_s = []
psi_r_d = []
psi_r = []

for idx, r in enumerate(rr):
    this_psi_r_s, this_psi_r_d = phi_r(r)
    psi_r_s.append(this_psi_r_s)
    psi_r_d.append(this_psi_r_d)

psi_r_s = np.array(psi_r_s)
psi_r_d = np.array(psi_r_d)
norm_psi = 1 / np.sqrt(np.sum(wr * rr**2 * ((psi_r_s) ** 2 + (psi_r_d) ** 2)))
psi_r_s = psi_r_s * norm_psi
psi_r_d = psi_r_d * norm_psi

psi_norm = np.sum(wr * rr**2 * ((psi_r_s) ** 2 + (psi_r_d) ** 2))
print(f"psi_norm = {psi_norm}")

r2 = np.sum(wr * rr**2 * (rr / 2) ** 2 * ((psi_r_s) ** 2 + (psi_r_d) ** 2))
print(f"r = {np.sqrt(r2)}")
Pd = np.sum(wr * rr**2 * (psi_r_d) ** 2)
print(f"Pd = {Pd}")
t5 = time.time()
profiler.add_timing("Get observables", t5 - t4)

################################################################################################################


utility.section_message("Timings")

profiler.print_timings()

################################################################################################################

utility.footer_message()
