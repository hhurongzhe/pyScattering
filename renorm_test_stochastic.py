import sys

sys.path.append("./lib")
import utility
import profiler
import time
import constants as const
import stochastic_srg as stochastic_srg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


params = {}
params["potential_type"] = "n3loemn500"
params["q_min"] = 1e-8
params["q_max"] = 5.0
params["q_number"] = 100
params["target_walker_number"] = 10000
params["random_sampling"] = False
params["loops"] = 1
params["d_tau"] = 1e-7
params["A"] = 1
params["xi"] = 0.1
params["zeta"] = 0.0025
params["initiator_approximation"] = False
params["initiator_threshold"] = 0.1
params["seed"] = 0


utility.header_message()

################################################################################################################

utility.section_message("Initialization")

s_target = 1e-1  # in fm^(4)
Lambda = pow(s_target, -1 / 4)  # in fm^(-1)
print("Lambda = ", Lambda, "fm^(-1)")
units_factor = const.MN**2 / const.hbarc**4
s_target *= units_factor


sSRG = stochastic_srg.sSRG(params)


sSRG.initialize_walkers()
sSRG.start(s_target)
mtx = sSRG.get_V()


plt.plot(figsize=(5, 5))
pp, p = np.meshgrid(sSRG.mesh_q, sSRG.mesh_q)
if mtx.min() >= 0:
    norm = TwoSlopeNorm(vmin=0, vmax=mtx.max())
elif mtx.max() <= 0:
    norm = TwoSlopeNorm(vmin=mtx.min(), vmax=0)
else:
    norm = TwoSlopeNorm(vmin=mtx.min(), vcenter=0, vmax=mtx.max())
# c = plt.imshow(mtx, cmap="RdBu_r", interpolation="bicubic", extent=(p.min(), p.max(), pp.min(), pp.max()), origin="lower", norm=norm)
c = plt.imshow(mtx, cmap="RdBu_r", interpolation="none", extent=(p.min(), p.max(), pp.min(), pp.max()), origin="lower", norm=norm)
plt.xlabel(r"$p$ (MeV)", fontsize=16)
plt.ylabel(r"$p'$ (MeV)", fontsize=16)
plt.title(r"$\lambda=$" + str(round(Lambda, 2)) + r"$\;\mathrm{fm}^{-1}$", fontsize=16)
plt.xticks([200, 400, 600, 800])
plt.yticks([200, 400, 600, 800])
plt.tick_params(labelsize=14)
cbar = plt.colorbar(c, orientation="vertical", pad=0.1, shrink=1)
cbar.set_label("$V(p',p)\,\mathrm{(MeV^{-2})}$", fontsize=16)
plt.tight_layout()
plt.savefig(f"srg-stochastic-single-channel-1s0-{params['potential_type']}.png", bbox_inches="tight", dpi=600)
plt.show()
