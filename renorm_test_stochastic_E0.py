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
params["flag"] = "3sd1"
params["coupled_channel"] = True
params["quantum_numbers"] = 1  # J
params["q_min"] = 1e-8
params["q_max"] = 5.0
params["q_number"] = 100
params["mesh_type"] = "linear"
params["target_walker_number"] = 4 * 1000
params["random_sampling"] = False
params["loops"] = 100
params["d_tau"] = 5e-13
params["A"] = 1
params["xi"] = 0.1
params["zeta"] = 0.0025
params["initiator_approximation"] = False
params["initiator_threshold"] = 0.1
params["seed"] = 0


utility.header_message()

################################################################################################################

utility.section_message("Initialization")

Lambda = 2.0
s_target = np.power(Lambda, -4)
print("Lambda = ", Lambda, "fm^(-1)")
print("s_target = ", s_target, "fm^(4)")
units_factor = const.MN**2 / const.hbarc**4
s_target *= units_factor


sSRG = stochastic_srg.sSRG(params)
print(f"E0 = {sSRG.E0}")

sSRG.initialize_walkers()
sSRG.start(s_target)
mean_tau, std_tau = sSRG.get_stat_array(sSRG.tau_loops_trace)
mean_E, std_E = sSRG.get_stat_array(sSRG.energy_loops_trace)

file_tau_name = f"result/srg-stoch-tau-Lambda{Lambda}.npy"
file_E0_mean_name = f"result/srg-stoch-E0-mean-{params['flag']}-{params['potential_type']}-Lambda{Lambda}-loop{params['loops']}-Nw{params['target_walker_number']}.npy"
file_E0_std_name = f"result/srg-stoch-E0-std-{params['flag']}-{params['potential_type']}-Lambda{Lambda}-loop{params['loops']}-Nw{params['target_walker_number']}.npy"
np.save(file_tau_name, mean_tau)
np.save(file_E0_mean_name, mean_E)
np.save(file_E0_std_name, std_E)


plt.plot(figsize=(5, 5))
plt.xlabel(r"$s$ [fm$^4$]", fontsize=16)
plt.ylabel(r"$E_0$ [MeV]", fontsize=16)
plt.tick_params(labelsize=14)
plt.plot(mean_tau, mean_E, c="C0", linestyle="-", linewidth=0.5)
plt.errorbar(mean_tau, mean_E, yerr=std_E, fmt="o", c="C0", markersize=3, linewidth=0.6, elinewidth=0.6, capsize=0, markeredgewidth=0, label=r"$E_0$ [MeV]")
plt.tight_layout()
plt.savefig(f"srg-stochastic-E0-coupled-channel-{params['flag']}-{params['potential_type']}.png", bbox_inches="tight", dpi=600)
plt.show()
