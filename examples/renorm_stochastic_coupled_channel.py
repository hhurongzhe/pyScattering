import sys

sys.path.append("./lib")
import os
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
params["q_max"] = 6.0
params["q_number"] = 100
params["mesh_type"] = "linear"
params["target_walker_number"] = 100
params["random_sampling"] = False
params["loops"] = 10
params["steps"] = 10000
params["A"] = 1
params["xi"] = 0.1
params["zeta"] = 0.0025
params["initiator_approximation"] = False
params["initiator_threshold"] = 0.1
params["seed"] = 0
params["phase_Tlabs"] = [1e-3] + [x / 10 for x in np.arange(1, 11, 0.1)] + [x for x in np.arange(2, 31, 1)] + [x for x in np.arange(40, 300, 10)]


utility.header_message()

################################################################################################################

utility.section_message("Initialization")


sSRG = stochastic_srg.sSRG(params)
sSRG.initialize_walkers()
sSRG.start()
mean_tau, std_tau = sSRG.get_stat_array(sSRG.tau_loops_trace)
mean_mtx, std_mtx = sSRG.get_stat_mtx()
mean_phase, std_phase = sSRG.get_phase_stat_array(sSRG.phase_loops_trace)
exact_mtx = sSRG.solve_exact_srg()
exact_phase = sSRG.compute_phase_shifts(exact_mtx)

os.makedirs("result", exist_ok=True)
file_tau_name = f"result/srg-stoch-tau-step{params['steps']}-Lambda{params['Lambda']}.npy"
file_mean_name = f"result/srg-stoch-mean-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
file_std_name = f"result/srg-stoch-std-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
file_exact_name = f"result/srg-exact-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
file_phase_raw_name = f"result/srg-stoch-phase-raw-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
file_phase_mean_name = f"result/srg-stoch-phase-mean-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
file_phase_std_name = f"result/srg-stoch-phase-std-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
file_phase_txt_name = f"result/srg-stoch-phase-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.txt"
file_exact_phase_name = f"result/srg-exact-phase-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.npy"
file_exact_phase_txt_name = f"result/srg-exact-phase-{params['flag']}-{params['potential_type']}-Lambda{params['Lambda']}-loop{params['loops']}-step{params['steps']}-Nw{params['target_walker_number']}.txt"
np.save(file_tau_name, mean_tau)
np.save(file_mean_name, mean_mtx)
np.save(file_std_name, std_mtx)
np.save(file_exact_name, exact_mtx)
np.save(file_phase_raw_name, np.array(sSRG.phase_loops_trace))
np.save(file_phase_mean_name, mean_phase)
np.save(file_phase_std_name, std_phase)
np.savetxt(
    file_phase_txt_name,
    np.column_stack((sSRG.phase_Tlabs, mean_phase, std_phase)),
    header="Tlab delta_minus_mean delta_plus_mean epsilon_mean delta_minus_std delta_plus_std epsilon_std",
)
np.save(file_exact_phase_name, exact_phase)
np.savetxt(
    file_exact_phase_txt_name,
    np.column_stack((sSRG.phase_Tlabs, exact_phase)),
    header="Tlab delta_minus_exact delta_plus_exact epsilon_exact",
)

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
cbar.set_label(r"$V(p',p)\,\mathrm{(MeV^{-2})}$", fontsize=16)
plt.tight_layout()
plt.savefig(f"srg-stochastic-coupled-channel-{params['flag']}-{params['potential_type']}.png", bbox_inches="tight", dpi=600)
plt.close()
