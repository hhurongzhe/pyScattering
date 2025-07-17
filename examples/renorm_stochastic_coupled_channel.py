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
params["Lambda"] = 2.0
params["flag"] = "3sd1"
params["coupled_channel"] = True
params["quantum_numbers"] = 1  # J
params["q_min"] = 1e-8
params["q_max"] = 5.0
params["q_number"] = 100
params["mesh_type"] = "linear"
params["target_walker_number"] = 1000
params["random_sampling"] = False
params["loops"] = 10
params["steps"] = 1000
params["A"] = 1
params["xi"] = 0.1
params["zeta"] = 0.0025
params["initiator_approximation"] = False
params["initiator_threshold"] = 0.1
params["seed"] = 0


utility.header_message()

################################################################################################################

utility.section_message("Initialization")


sSRG = stochastic_srg.sSRG(params)
sSRG.initialize_walkers()
sSRG.start()
mean_tau, std_tau = sSRG.get_stat_array(sSRG.tau_loops_trace)
mean_mtx, std_mtx = sSRG.get_stat_mtx()

file_tau_name = f"result/srg-stoch-tau-step{params['steps']}-Lambda{params['Lambda']}.npy"
file_mean_name = f"result/srg-stoch-mean-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
file_std_name = f"result/srg-stoch-std-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
np.save(file_tau_name, mean_tau)
np.save(file_mean_name, mean_mtx)
np.save(file_std_name, std_mtx)

plt.plot(figsize=(5, 5))
pp, p = np.meshgrid(sSRG.mesh_q, sSRG.mesh_q)
if mean_mtx.min() >= 0:
    norm = TwoSlopeNorm(vmin=0, vmax=mean_mtx.max())
elif mean_mtx.max() <= 0:
    norm = TwoSlopeNorm(vmin=mean_mtx.min(), vmax=0)
else:
    norm = TwoSlopeNorm(vmin=mean_mtx.min(), vcenter=0, vmax=mean_mtx.max())
# c = plt.imshow(mtx, cmap="RdBu_r", interpolation="bicubic", extent=(p.min(), p.max(), pp.min(), pp.max()), origin="lower", norm=norm)
c = plt.imshow(mean_mtx, cmap="RdBu_r", interpolation="none", extent=(p.min(), p.max(), pp.min(), pp.max()), origin="lower", norm=norm)
plt.xlabel(r"$p$ (MeV)", fontsize=16)
plt.ylabel(r"$p'$ (MeV)", fontsize=16)
plt.title(r"$\lambda=$" + str(round(params["Lambda"], 2)) + r"$\;\mathrm{fm}^{-1}$", fontsize=16)
plt.xticks([200, 400, 600, 800])
plt.yticks([200, 400, 600, 800])
plt.tick_params(labelsize=14)
cbar = plt.colorbar(c, orientation="vertical", pad=0.1, shrink=1)
cbar.set_label("$V(p',p)\,\mathrm{(MeV^{-2})}$", fontsize=16)
plt.tight_layout()
plt.savefig(f"srg-stochastic-coupled-channel-{params['flag']}-{params['potential_type']}.png", bbox_inches="tight", dpi=600)
plt.show()
