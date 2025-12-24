import sys

sys.path.append("./lib")

import numpy as np
import nn_studio as nn_studio
import chiral_potential as chiral_potential
import matplotlib.pyplot as plt

# initialize an object for computing T-matrices, phase shifts,
nn = nn_studio.nn_studio(jmin=0, jmax=1, tz=0, Np=128, mesh_type="gauleg_finite")

potential_type = "n3loem"

# initialize an object for the chiral interaction
potential = chiral_potential.two_nucleon_potential(potential_type)

# give the potential to the nn-analyzer
nn.V = potential

is_coupled, l, s, j, t, tz = False, 0, 0, 0, 0, 0
potential_matrix = nn.setup_Vmtx(is_coupled, (l, s, j, t, tz))

pp, p = np.meshgrid(nn.pmesh, nn.pmesh)


z_min, z_max = -np.abs(potential_matrix).max(), np.abs(potential_matrix).max()
fig, ax = plt.subplots(figsize=(4, 3))
c = ax.pcolormesh(p, pp, potential_matrix, cmap="RdBu", vmin=z_min, vmax=z_max)
fig.colorbar(c, ax=ax)
ax.set_xlabel(r"$p$ (MeV)")
ax.set_ylabel(r"$p'$ (MeV)")
plt.savefig(f"plot_potential_nonlocal_{potential_type}.png", bbox_inches="tight", dpi=600)
plt.show()
