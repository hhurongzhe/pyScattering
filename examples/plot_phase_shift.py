import numpy as np
import matplotlib.pyplot as plt
from lib.ini_plot import *


def read_phase_shift_data(file_path, name):
    with open(file_path, "r") as f:
        headers = f.readline().strip().split()
    if name not in headers:
        raise ValueError(f"列名 '{name}' 不存在。可选列名: {headers}")
    col_idx = headers.index(name)
    data = np.loadtxt(file_path, skiprows=1)
    if name == "3D3":  # 单独处理一下3D3的相移
        return np.where(data[:, col_idx] > 90, 180 - data[:, col_idx], data[:, col_idx])
    return data[:, col_idx]


chan = "np"
int_names = ["av18"]
int_labels = [r"$\mathrm{av18}$"]
int_colors = ["C0", "C1", "C2", "C3", "C4", "C9"]
pw_labels = [r"^1S_0", r"^3P_0", r"\epsilon_1", r"^1D_2", r"^1P_1", r"^3P_1", r"^3D_2", r"^3P_2", r"^3S_1", r"^3D_1", r"\epsilon_2", r"^3D_3"]
pw_names = ["1S0", "3P0", "E1", "1D2", "1P1", "3P1", "3D2", "3P2", "3S1", "3D1", "E2", "3D3"]

ini_plot()
fig, axs = plt.subplots(4, 3, figsize=(4, 5))
fig.subplots_adjust(hspace=0.2, wspace=0.25)
fs = 6
for i, pw in enumerate(pw_labels):
    ax = axs.flatten()[i]
    x_pos, y_pos = 0.5, 0.8
    if pw in [r"^3P_0"]:
        x_pos, y_pos = 0.2, 0.2
    if pw in [r"^1D_2", r"^3D_2", r"^3P_2"]:
        x_pos, y_pos = 0.15, 0.75
    ax.text(x_pos, y_pos, rf"${pw}$", fontsize=fs + 2, transform=ax.transAxes, ha="center")
    # ax.set_xlabel(r"$T_\mathrm{lab}$ [MeV]", fontsize=fs)
    # ax.set_ylabel(r"$\delta$ [deg]", fontsize=fs)
    ax.minorticks_on()
    ax.tick_params(direction="in", width=0.5, length=2, top=True, bottom=True, left=True, right=True, axis="both", labelsize=5)
    ax.tick_params(which="minor", direction="in", width=0.3, length=1, top=True, bottom=True, left=True, right=True, axis="both")
    ax.tick_params(which="minor", length=2, width=0.3, axis="x")
    # plot nijmegen and granada data
    phase_nij = read_phase_shift_data(f"lib/data_phase/phase_shifts_nijmegen_{chan}.txt", pw_names[i])
    Tlabs_nij = read_phase_shift_data(f"lib/data_phase/phase_shifts_nijmegen_{chan}.txt", "Tlab")
    phase_granada = read_phase_shift_data(f"lib/data_phase/phase_shifts_granada_{chan}.txt", pw_names[i])
    error_granada = read_phase_shift_data(f"lib/data_phase/phase_shifts_granada_{chan}.txt", pw_names[i] + "_err")
    Tlabs_granada = read_phase_shift_data(f"lib/data_phase/phase_shifts_granada_{chan}.txt", "Tlab")
    # plot calculated phase shifts
    for j, int_name in enumerate(int_names):
        file_path = f"result/phase_shifts_{int_name}_Jmax3_{chan}.txt"
        phase_shifts = read_phase_shift_data(file_path, pw_names[i])
        Tlabs = read_phase_shift_data(file_path, "Tlab")
        ax.plot(Tlabs, phase_shifts, color=int_colors[j], linewidth=0.7)
    ax.plot(Tlabs_nij, phase_nij, "D", markersize=2, markeredgewidth=0.6, markeredgecolor="black", markerfacecolor="none", zorder=2)
    ax.errorbar(Tlabs_granada, phase_granada, yerr=error_granada, fmt="o", c="black", capsize=0, elinewidth=0.6, capthick=0.6, markersize=2.5, markeredgewidth=0.6, markeredgecolor="black", markerfacecolor="none", zorder=2)

# 设置公用图例
ax = axs.flatten()[-1]
for j, int_name in enumerate(int_names):
    ax.plot([], [], label=int_labels[j], color=int_colors[j], linewidth=0.7)
ax.plot([], [], "D", markersize=2, markeredgewidth=0.6, markeredgecolor="black", markerfacecolor="none", label="Nijmegen")
ax.errorbar([], [], yerr=[], fmt="o", c="black", capsize=0, elinewidth=0.6, capthick=0.6, markersize=2.5, markeredgewidth=0.6, markeredgecolor="black", markerfacecolor="none", label=r"Granada")
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.51, 0.88), fontsize=fs + 1)

fig.supxlabel(r"$T_\mathrm{lab}$ [MeV]", fontsize=fs + 2, y=0.05)
fig.supylabel(r"$\delta$ [deg]", fontsize=fs + 2, x=0.05)

for spine in ["bottom", "top", "left", "right"]:
    for ax in axs.flatten():
        ax.spines[spine].set_linewidth(0.5)
# plt.savefig("fig_phase_shift.pdf", bbox_inches="tight", transparent=False, dpi=1200)
plt.show()
