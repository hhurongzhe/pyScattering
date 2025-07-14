import sys

sys.path.append("./lib")

import numpy as np
import chiral_potential as chiral_potential
import matplotlib.pyplot as plt


# initialize an object for the chiral interaction
potential1 = chiral_potential.two_nucleon_potential("idaholocal12")

r = np.linspace(1e-6, 4, 200)

ll, l, s, j, tz = 1, 1, 1, 0, 0
mtx1 = [potential1.potential_local(ll, l, s, j, tz, rr) for rr in r]

fig = plt.figure(figsize=(6, 6), dpi=160)
plt.xlim(0, 4)
# plt.ylim(-200, 200)
plt.xlabel(r"$r (fm)$")
plt.ylabel(r"$V(r) (\mathrm{MeV})$")
plt.plot(r, mtx1, color="C2", linestyle="-", linewidth=1.2)
plt.savefig("plot_potential_local.png", bbox_inches="tight", dpi=600)
plt.show()
