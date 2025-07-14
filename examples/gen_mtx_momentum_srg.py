import sys

sys.path.append("./lib")
import os
import struct
import utility
import profiler
import time
import chiral_potential as chiral_potential
import numpy as np
import constants as const
from scipy.integrate import solve_ivp

mp = 938.2720
mn = 939.5654
mN = 2 * (mp * mn) / (mp + mn)

potential_type = "n3loem"
potential = chiral_potential.two_nucleon_potential(potential_type)
print(f"potential type : {potential_type}\n")

S_eval_points = np.array([0.0952599])
Ns = len(S_eval_points)

# from fm^(4) to MeV^(-2).
units_factor = mN**2 / const.hbarc**4
S_eval_points = units_factor * S_eval_points

# Total number of integration points Nq~100 would be enough for most cases.
# For Chiral EFT interactions, use q_max 6-8 fm^(-1).
# For standard OBE interactions such as AV18, use q_max 20-30 fm^(-1).
Nq = 100
q_max = 8
q_min = 1e-14


def gauss_legendre_line_mesh(a, b):
    x, w = np.polynomial.legendre.leggauss(Nq)
    # Translate x values from the interval [-1, 1] to [a, b]
    t = 0.5 * (x + 1) * (b - a) + a
    u = w * 0.5 * (b - a)

    return t, u


mesh_q, weight_q = gauss_legendre_line_mesh(q_min * const.hbarc, q_max * const.hbarc)


def gen_mtx_srg_coupled(J: int, Tz: int):

    def potential_mm_initial(pp: float, p: float):
        ll, l, j, s, tz = J - 1, J - 1, J, 1, Tz
        return potential.potential(ll, l, pp, p, j, s, tz)

    def potential_mp_initial(pp: float, p: float):
        ll, l, j, s, tz = J - 1, J + 1, J, 1, Tz
        return potential.potential(ll, l, pp, p, j, s, tz)

    def potential_pm_initial(pp: float, p: float):
        ll, l, j, s, tz = J + 1, J - 1, J, 1, Tz
        return potential.potential(ll, l, pp, p, j, s, tz)

    def potential_pp_initial(pp: float, p: float):
        ll, l, j, s, tz = J + 1, J + 1, J, 1, Tz
        return potential.potential(ll, l, pp, p, j, s, tz)

    V0_mm = np.zeros((Nq, Nq))
    V0_mp = np.zeros((Nq, Nq))
    V0_pm = np.zeros((Nq, Nq))
    V0_pp = np.zeros((Nq, Nq))
    Tkin = np.zeros((Nq, Nq))
    for idx_p, p in enumerate(mesh_q):
        Tkin[idx_p, idx_p] = p**2 / mN
        for idx_pp, pp in enumerate(mesh_q):
            w_p = weight_q[idx_p]
            w_pp = weight_q[idx_pp]
            V0_mm[idx_pp, idx_p] = np.sqrt(w_pp) * pp * potential_mm_initial(pp, p) * np.sqrt(w_p) * p
            V0_mp[idx_pp, idx_p] = np.sqrt(w_pp) * pp * potential_mp_initial(pp, p) * np.sqrt(w_p) * p
            V0_pm[idx_pp, idx_p] = np.sqrt(w_pp) * pp * potential_pm_initial(pp, p) * np.sqrt(w_p) * p
            V0_pp[idx_pp, idx_p] = np.sqrt(w_pp) * pp * potential_pp_initial(pp, p) * np.sqrt(w_p) * p
    Tkin_coupled = np.block([[Tkin, np.zeros((Nq, Nq))], [np.zeros((Nq, Nq)), Tkin]])
    V0_coupled = np.block([[V0_mm, V0_mp], [V0_pm, V0_pp]])

    print("initialize V0 complete.")

    # commutator, arguments A and B matrices
    def commutator(A, B):
        return np.matmul(A, B) - np.matmul(B, A)

    # commutator for coupled channels including reshaping to [Nrows,Nrows] matrices, arguments A and B are 1d arrays
    def commutator_reshape_coupled(A: np.ndarray, B: np.ndarray, Nrows: int):
        A_reshape = np.reshape(A, (2 * Nrows, 2 * Nrows))
        B_reshape = np.reshape(B, (2 * Nrows, 2 * Nrows))
        return np.matmul(A_reshape, B_reshape) - np.matmul(B_reshape, A_reshape)

    def dVds_coupled(s: float, Vcoupled: np.ndarray, Nrows: int, Tkincoupled: np.ndarray):
        Vcoupled_reshape = np.reshape(Vcoupled, (2 * Nrows, 2 * Nrows))
        ret = np.reshape(
            commutator_reshape_coupled(
                commutator_reshape_coupled(Tkincoupled, Vcoupled_reshape, Nrows),
                Tkincoupled + Vcoupled_reshape,
                Nrows,
            ),
            Vcoupled.size,
        )
        return ret

    solution = solve_ivp(
        dVds_coupled,
        [0, 1.1 * S_eval_points.max()],
        V0_coupled.ravel(),
        method="RK23",
        t_eval=S_eval_points,
        args=(Nq, Tkin_coupled),
        rtol=1e-8,
        atol=1e-8,
    )
    if not solution.success:
        error_message = "The ODE solver failed with error message:\n" + solution.message
        sys.exit(error_message)
    print("Solving SRG flow complete.")

    idx_s, s_plot = 0, S_eval_points[0]
    S_to_plot = s_plot / units_factor  # in fm^(4).
    Lambda = str(round(pow(S_to_plot, -1 / 4), 1))

    potential_plot = solution.y[:, idx_s].reshape((2 * Nq, 2 * Nq))

    potential_plot_without_factor = np.zeros_like(V0_coupled)
    for idxpp, pp in enumerate(mesh_q):
        for idxp, p in enumerate(mesh_q):
            w_p = weight_q[idxp]
            w_pp = weight_q[idxpp]
            potential_plot_without_factor[idxpp][idxp] = potential_plot[idxpp][idxp] / (np.sqrt(w_pp * w_p) * pp * p)
            potential_plot_without_factor[idxpp + Nq][idxp] = potential_plot[idxpp + Nq][idxp] / (np.sqrt(w_pp * w_p) * pp * p)
            potential_plot_without_factor[idxpp][idxp + Nq] = potential_plot[idxpp][idxp + Nq] / (np.sqrt(w_pp * w_p) * pp * p)
            potential_plot_without_factor[idxpp + Nq][idxp + Nq] = potential_plot[idxpp + Nq][idxp + Nq] / (np.sqrt(w_pp * w_p) * pp * p)
    return potential_plot_without_factor


def gen_mtx_srg_single(L: int, S: int, J: int, Tz: int):
    def potential_single_initial(pp: float, p: float):
        ll, l, j, s, tz = L, L, J, S, Tz
        return potential.potential(ll, l, pp, p, j, s, tz)

    V0 = np.zeros((Nq, Nq))
    Tkin = np.zeros((Nq, Nq))
    for idx_p, p in enumerate(mesh_q):
        Tkin[idx_p, idx_p] = p**2 / mN
        for idx_pp, pp in enumerate(mesh_q):
            w_p = weight_q[idx_p]
            w_pp = weight_q[idx_pp]
            V0[idx_pp, idx_p] = np.sqrt(w_pp) * pp * potential_single_initial(pp, p) * np.sqrt(w_p) * p
    print("initialize V0 complete.")

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
            commutator_reshape(commutator_reshape(Tkin, V_reshape, Nrows), Tkin + V_reshape, Nrows),
            V.size,
        )
        return ret

    solution = solve_ivp(
        dVds,
        [0, 1.1 * S_eval_points.max()],
        V0.ravel(),
        method="RK23",
        t_eval=S_eval_points,
        args=(Nq, Tkin),
        rtol=1e-8,
        atol=1e-8,
    )
    if not solution.success:
        error_message = "The ODE solver failed with error message:\n" + solution.message
        sys.exit(error_message)
    print("Solving SRG flow complete.")

    idx_s, s_plot = 0, S_eval_points[0]
    S_to_plot = s_plot / units_factor  # in fm^(4).
    Lambda = str(round(pow(S_to_plot, -1 / 4), 1))

    potential_plot = solution.y[:, idx_s].reshape((Nq, Nq))

    potential_plot_without_factor = np.zeros_like(V0)
    for idxpp, pp in enumerate(mesh_q):
        for idxp, p in enumerate(mesh_q):
            w_p = weight_q[idxp]
            w_pp = weight_q[idxpp]
            potential_plot_without_factor[idxpp][idxp] = potential_plot[idxpp][idxp] / (np.sqrt(w_pp * w_p) * pp * p)
    return potential_plot_without_factor


def gen_mtx_channels(N, Jmax):
    assert Jmax >= 0, "Jmax must be nonegative"

    # pp channels first
    tz = -1
    channels_pp = []
    channels_pp.append([0, 1, 0, tz, N])
    channels_pp.append([0, -1, 1, tz, N])
    if Jmax > 0:
        for j in range(1, Jmax + 1, 1):
            if j % 2 == 1:
                channels_pp.append([j, -1, 1, tz, N])
            else:
                channels_pp.append([j, 1, 0, tz, N])
                channels_pp.append([j, -1, 1, tz, 2 * N])

    # pn channels second
    tz = 0
    channels_pn = []
    channels_pn.append([0, 1, 0, tz, N])
    channels_pn.append([0, -1, 1, tz, N])
    if Jmax > 0:
        for j in range(1, Jmax + 1, 1):
            if j % 2 == 1:
                channels_pn.append([j, 1, 1, tz, 2 * N])
                channels_pn.append([j, -1, 0, tz, N])
                channels_pn.append([j, -1, 1, tz, N])
            else:
                channels_pn.append([j, 1, 0, tz, N])
                channels_pn.append([j, 1, 1, tz, N])
                channels_pn.append([j, -1, 1, tz, 2 * N])

    # nn channels last
    tz = 1
    channels_nn = []
    channels_nn.append([0, 1, 0, tz, N])
    channels_nn.append([0, -1, 1, tz, N])
    if Jmax > 0:
        for j in range(1, Jmax + 1, 1):
            if j % 2 == 1:
                channels_nn.append([j, -1, 1, tz, N])
            else:
                channels_nn.append([j, 1, 0, tz, N])
                channels_nn.append([j, -1, 1, tz, 2 * N])
    return channels_pp + channels_pn + channels_nn


def convert_dat_to_binary(input_file, output_file):
    with open(input_file, "r") as fp_in, open(output_file, "wb") as fp_out:
        lines = fp_in.readlines()

        # Parse and write NMesh, Jmax, and NChan to binary file
        NMesh = int(lines[1])
        Jmax = int(lines[3])
        NChan = int(lines[5])
        fp_out.write(struct.pack("<i", NMesh))
        fp_out.write(struct.pack("<i", Jmax))
        fp_out.write(struct.pack("<i", NChan))

        # Parse and write Momentum Mesh Points to binary file
        start_index = 7
        for i in range(NMesh):
            fp_out.write(struct.pack("<d", float(lines[start_index + i])))

        # Parse and write Momentum Mesh Weights to binary file
        start_index = start_index + NMesh + 1
        for i in range(NMesh):
            fp_out.write(struct.pack("<d", float(lines[start_index + i])))

        # Parse and write J, Prty, S, Tz, Ndim, and V values to binary file
        start_index = 2 * NMesh + 9
        for _ in range(NChan):
            J = int(lines[start_index])
            Prty = int(lines[start_index + 2])
            S = int(lines[start_index + 4])
            Tz = int(lines[start_index + 6])
            Ndim = int(lines[start_index + 8])
            fp_out.write(struct.pack("<i", J))
            fp_out.write(struct.pack("<i", Prty))
            fp_out.write(struct.pack("<i", S))
            fp_out.write(struct.pack("<i", Tz))
            fp_out.write(struct.pack("<i", Ndim))

            V_lines = lines[start_index + 10 : start_index + 10 + Ndim]
            V_values = [float(value) for line in V_lines for value in line.split()]
            for value in V_values:
                fp_out.write(struct.pack("<d", value))

            start_index += 11 + Ndim


# Writter of interaction matrix elements in momentum space, where [kmax] = fm^(-1).
# The resulting interaction mtx are in [MeV fm^(-3)].
def write_mtx_momentum(Jmax=8):
    t1 = time.time()
    # this value of hw is in line with Miyagi's code,
    # used for unit transformation.
    hc = 197.32705
    hccubic = hc**3
    dir = "input_nn_files"  # files are generated in this dir.
    if not os.path.exists(dir):
        os.makedirs(dir)
    file_name = "./" + dir + "/" + f"{potential_type}_kmin{q_min}_kmax{q_max}_N{Nq}_Jmax{Jmax}.dat"
    channels = gen_mtx_channels(Nq, Jmax)
    MeshPoints, MeshWeights = gauss_legendre_line_mesh(q_min, q_max)
    with open(file_name, "w") as fp:
        fp.write(f"NMesh:\n{Nq}\n")
        fp.write(f"Jmax:\n{Jmax}\n")
        fp.write(f"NChan:\n{len(channels)}\n")
        fp.write("Momentum Mesh Points:\n")
        np.savetxt(fp, MeshPoints, fmt="%.17f")
        fp.write("Momentum Mesh Weights:\n")
        np.savetxt(fp, MeshWeights, fmt="%.17f")
        MeshPoints = MeshPoints * hc
        MeshWeights = MeshWeights * hc
        for chan in channels:
            J = chan[0]
            Prty = chan[1]
            S = chan[2]
            Tz = chan[3]
            Ndim = chan[4]
            coup = False
            if Ndim != Nq:
                coup = True
            fp.write(f"J:\n{J}\n")
            fp.write(f"Prty:\n{Prty}\n")
            fp.write(f"S:\n{S}\n")
            fp.write(f"Tz:\n{Tz}\n")
            fp.write(f"Ndim:\n{Ndim}\n")
            fp.write("V:\n")
            if J == 0:
                V = gen_mtx_srg_single(S, S, J, Tz) * hccubic
                np.savetxt(fp, V, fmt="%.17f")
            else:
                if not coup:
                    V = gen_mtx_srg_single(J, S, J, Tz) * hccubic
                    np.savetxt(fp, V, fmt="%.17f")
                else:
                    V = gen_mtx_srg_coupled(J, Tz) * hccubic
                    np.savetxt(fp, V, fmt="%.17f")
    binary_file = "./" + dir + "/" + f"{potential_type}_kmin{q_min}_kmax{q_max}_N{Nq}_Jmax{Jmax}.bin"
    t2 = time.time()
    profiler.add_timing("Writing .dat", t2 - t1)
    convert_dat_to_binary(file_name, binary_file)
    t3 = time.time()
    profiler.add_timing("Writing .bin", t3 - t2)
    print("files are generated to :\n", file_name, "\n", binary_file)


write_mtx_momentum(Jmax=8)
