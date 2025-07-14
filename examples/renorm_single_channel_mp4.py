## a srg solver in momentum space,
## written by hu rongzhe.

import sys

sys.path.append("./lib")
import utility
import profiler
import time
import chiral_potential as chiral_potential
import numpy as np
import matplotlib.pyplot as plt
import constants as const
from scipy.integrate import solve_ivp


utility.header_message()

################################################################################################################

utility.section_message("Initialization")

mp = 938.2720
mn = 939.5654
mN = 2 * (mp * mn) / (mp + mn)

potential_type = "n3loem"
potential = chiral_potential.two_nucleon_potential(potential_type)
print(f"potential type : {potential_type}\n")

# SRG flow parameters which you want to get, in fm^(4).
Ns = 50
S_eval_points = np.logspace(-4, -0.6, Ns)


# from fm^(4) to MeV^(-2).
units_factor = mN**2 / const.hbarc**4
S_eval_points = units_factor * S_eval_points

# Total number of integration points Nq~100 would be enough for most cases.
# For Chiral EFT interactions, use q_max 6-8 fm^(-1).
# For standard OBE interactions such as AV18, use q_max 20-30 fm^(-1).
Nq = 200
q_max = 5.0


# we use linear mesh points.
q_min = 1e-8
Delta_q = (q_max - q_min) / Nq * const.hbarc
mesh_q = const.hbarc * np.linspace(q_min, q_max, Nq)
weight_q = np.full(Nq, Delta_q)


# np interaction in 1S0 channel.
def potential_1s0_initial(pp: float, p: float):
    ll, l, j, s, tz = 0, 0, 0, 0, 0
    return potential.potential(ll, l, pp, p, j, s, tz)


t1 = time.time()
V0 = np.zeros((Nq, Nq))
Tkin = np.zeros((Nq, Nq))
for idx_p, p in enumerate(mesh_q):
    Tkin[idx_p, idx_p] = p**2 / mN
    for idx_pp, pp in enumerate(mesh_q):
        w_p = weight_q[idx_p]
        w_pp = weight_q[idx_pp]
        V0[idx_pp, idx_p] = (
            np.sqrt(w_pp) * pp * potential_1s0_initial(pp, p) * np.sqrt(w_p) * p
        )

t2 = time.time()
profiler.add_timing("Initialize V0", t2 - t1)
print("initialize V0 complete.")

################################################################################################################


################################################################################################################
utility.section_message("Solving SRG Flow Equations With Scipy")


# commutator, arguments A and B matrices
def commutator(A, B):
    return np.matmul(A, B) - np.matmul(B, A)


# commutator including reshaping to [Nrows,Nrows] matrices, arguments A and B are 1d arrays
def commutator_reshape(A, B, Nrows):
    A_reshape = np.reshape(A, (Nrows, Nrows))
    B_reshape = np.reshape(B, (Nrows, Nrows))
    return np.matmul(A_reshape, B_reshape) - np.matmul(B_reshape, A_reshape)


# right hand side of SRG flow equation, implementation of double commutators for uncoupled and coupled cases
def dVds(s: float, V, Nrows, Tkin):
    V_reshape = np.reshape(V, (Nrows, Nrows))
    ret = np.reshape(
        commutator_reshape(
            commutator_reshape(Tkin, V_reshape, Nrows), Tkin + V_reshape, Nrows
        ),
        V.size,
    )
    return ret


t3 = time.time()
solution = solve_ivp(
    dVds,
    [0, 1.1 * S_eval_points.max()],
    V0.ravel(),
    method="RK23",
    t_eval=S_eval_points,
    args=(Nq, Tkin),
    rtol=1e-5,
    atol=1e-5,
)
if not solution.success:
    error_message = "The ODE solver failed with error message:\n" + solution.message
    sys.exit(error_message)
t4 = time.time()
profiler.add_timing("SRG Flow", t4 - t3)
print("Solving SRG flow complete.")

################################################################################################################


################################################################################################################
utility.section_message("Plot")

import os
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import TwoSlopeNorm

# 设置FFmpeg的路径
os.environ["PATH"] += os.pathsep + r"E:\ffmpeg\bin"

fig, ax = plt.subplots(figsize=(10, 10), dpi=300)  # 增加图像尺寸和DPI
c = ax.pcolormesh(mesh_q, mesh_q, np.abs(V0), cmap="afmhot")
ax.set_xlabel(r"$p$ (MeV)")
ax.set_ylabel(r"$p'$ (MeV)")


# Function to update the plot for each frame
def update(frame):
    ax.clear()
    s_plot = S_eval_points[frame]
    S_to_plot = s_plot / units_factor  # in fm^(4).
    Lambda = str(round(pow(S_to_plot, -1 / 4), 1))

    potential_plot = solution.y[:, frame].reshape((Nq, Nq))

    potential_plot_without_factor = np.zeros_like(V0)
    for idxpp, pp in enumerate(mesh_q):
        for idxp, p in enumerate(mesh_q):
            w_p = weight_q[idxp]
            w_pp = weight_q[idxpp]
            potential_plot_without_factor[idxpp][idxp] = potential_plot[idxpp][idxp] / (
                np.sqrt(w_pp * w_p) * pp * p
            )

    mtx = potential_plot_without_factor
    norm = TwoSlopeNorm(vmin=mtx.min(), vcenter=0, vmax=mtx.max())
    c = ax.pcolormesh(
        mesh_q,
        mesh_q,
        mtx,
        cmap="RdBu_r",
        norm=norm,
    )
    ax.set_xlabel(r"$p$ (MeV)")
    ax.set_ylabel(r"$p'$ (MeV)")
    ax.set_title(r"$\lambda=$" + Lambda + r"$\;\mathrm{fm}^{-1}$")


# Create the animation
animation = FuncAnimation(
    fig, update, frames=len(S_eval_points), interval=200, repeat=True
)

# Save the animation as high-quality MP4
writer = FFMpegWriter(
    fps=16, metadata=dict(artist="Me"), bitrate=4000, extra_args=["-crf", "18"]
)  # 提高bitrate
animation.save(
    f"srg-single-channel-1s0-{potential_type}.mp4", writer=writer, dpi=300  # 提高DPI
)

# Show the animation
plt.show()


################################################################################################################

################################################################################################################

utility.section_message("Timings")

profiler.print_timings()

################################################################################################################

utility.footer_message()
