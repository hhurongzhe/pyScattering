import sys

sys.path.append("./lib")
import utility
import profiler
import time
import chiral_potential as chiral_potential
import basic_math as bm
import numpy as np
import math
import WignerSymbol as ws
from functools import lru_cache
from tqdm import tqdm

utility.header_message()

################################################################################################################

utility.section_message("Setting Parameters")

Mp = 938.27231
Mn = 939.56563
hc = 197.32705

potential_type = "n2loopt"
include_coulomb = False
hw = 16
emax = 2
e2max = 4
JMAX = 8
verbose = True
p_min, p_max, p_meshnumber = 0.0, 8.0 * hc, 100  # MeV; cutoff = 8 fm^-1

print(f"    potential type    : {potential_type}")
print(f"    include Coulomb   : {include_coulomb}")
print(f"    hw                : {hw}")
print(f"    emax              : {emax}")
print(f"    e2max             : {e2max}")
print(f"    Jmax              : {JMAX}")
print(f"    p min             : {p_min:.1f} MeV")
print(f"    p max             : {p_max:.5f} MeV")
print(f"    p mesh numbers    : {p_meshnumber}")

################################################################################################################


################################################################################################################

utility.section_message("Initialization")

t1 = time.time()

mu = Mp * Mn / (Mp + Mn)
alpha = np.sqrt(mu * hw)
potential = chiral_potential.two_nucleon_potential(potential_type)
print(f"potential type : {potential_type}\n")
p_points, p_weights = bm.gauss_legendre_line_mesh(p_min, p_max, p_meshnumber)
ws.init(24, "Jmax", 9)

t2 = time.time()
profiler.add_timing("Storing 6,9-j couplings", t2 - t1)
################################################################################################################


################################################################################################################
# defines necessary functions.


# Momentum-space LSJ block in MeV^-2; rows are p', columns are p.
@lru_cache(maxsize=None)
def v_LSJ_CM(ll: int, l: int, s: int, j: int, tz: int):
    return np.array([[potential.potential(ll, l, pp, p, j, s, tz) for p in p_points] for pp in p_points])


# Momentum HO radial function, normalized with p^2 dp, in MeV^(-3/2).
# (-1)^n matches the radial Fourier phase used with these LSJ potentials.
@lru_cache(maxsize=None)
def HO_wf_momentum(n: int, l: int, a: float, p: float):
    x = p / a
    norm = np.sqrt(2 * math.factorial(n) / math.gamma(n + l + 1.5)) / a**1.5
    return (-1) ** n * norm * np.exp(-0.5 * x * x) * x**l * bm.laguerre_poly(n, l + 0.5, x * x)


@lru_cache(maxsize=None)
def HO_overlap(n: int, l: int, a: float):
    return p_weights * p_points**2 * np.array([HO_wf_momentum(n, l, a, p) for p in p_points])


@lru_cache(maxsize=None)
def is_triangle(j1: float, j2: float, j3: float):
    precision = 1e-6
    boo1 = np.abs((j1 + j2 + j3) - int(j1 + j2 + j3)) < precision
    boo2 = int(np.abs(j1 - j2) - j3) <= 0
    boo3 = int(np.abs(j1 + j2) - j3) >= 0
    boo = boo1 and boo2 and boo3
    return boo


@lru_cache(maxsize=None)
def cg(j1, m1, j2, m2, j, m):
    dj1 = int(2.0 * j1)
    dj2 = int(2.0 * j2)
    dj = int(2.0 * j)
    dm1 = int(2.0 * m1)
    dm2 = int(2.0 * m2)
    dm = int(2.0 * m)
    result = ws.CG(dj1, dj2, dj, dm1, dm2, dm)
    return result


@lru_cache(maxsize=None)
def f6j(j1, j2, j3, j4, j5, j6):
    dj1 = int(2.0 * j1)
    dj2 = int(2.0 * j2)
    dj3 = int(2.0 * j3)
    dj4 = int(2.0 * j4)
    dj5 = int(2.0 * j5)
    dj6 = int(2.0 * j6)
    result = ws.f6j(dj1, dj2, dj3, dj4, dj5, dj6)
    return result


@lru_cache(maxsize=None)
def f9j(j1, j2, j3, j4, j5, j6, j7, j8, j9):
    dj1 = int(2.0 * j1)
    dj2 = int(2.0 * j2)
    dj3 = int(2.0 * j3)
    dj4 = int(2.0 * j4)
    dj5 = int(2.0 * j5)
    dj6 = int(2.0 * j6)
    dj7 = int(2.0 * j7)
    dj8 = int(2.0 * j8)
    dj9 = int(2.0 * j9)
    result = ws.f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)
    return result


@lru_cache(maxsize=None)
def moshinsky(
    N: int,
    L: int,
    n: int,
    l: int,
    n1: int,
    l1: int,
    n2: int,
    l2: int,
    Lambda: int,
    tan_beta: float,
):
    t1 = time.time()
    result = ws.Moshinsky(N, L, n, l, n1, l1, n2, l2, Lambda, tan_beta)
    t2 = time.time()
    profiler.add_timing("Cal Moshinsky", t2 - t1)
    return result


@lru_cache(maxsize=None)
def coulomb_radial_overlap(n: int, l: int, frequency: float):
    # pp reduced mass should be used here.
    r, w = bm.gauss_legendre_line_mesh(0.0, 40.0, 500)
    nu = np.sqrt(0.5 * Mp * frequency) / hc
    x = nu * r
    norm = np.sqrt(2 * nu * math.factorial(n) / math.gamma(n + l + 1.5))
    laguerre = np.array([bm.laguerre_poly(n, l + 0.5, xi * xi) for xi in x])
    u = norm * x ** (l + 1) * np.exp(-0.5 * x * x) * laguerre
    # u(r)=r R(r); coordinate-space HO functions have no (-1)^n phase.
    return u * np.sqrt(w * hc / (137.035999 * r))


@lru_cache(maxsize=None)
def coulomb_harmonic_CM(nn: int, n: int, l: int, frequency: float):
    return float(np.sum(coulomb_radial_overlap(nn, l, frequency) * coulomb_radial_overlap(n, l, frequency)))


# potential in C.M. frame under HO basis.
@lru_cache(maxsize=None)
def v_harmonic_CM(nn: int, n: int, ll: int, l: int, s: float, j: float, t: int, a: float):
    t1 = time.time()
    bra = HO_overlap(nn, ll, a)
    ket = HO_overlap(n, l, a)
    temp = float(bra @ v_LSJ_CM(ll, l, s, j, t) @ ket)
    if include_coulomb and t == -1 and ll == l and (l + s) % 2 == 0:
        temp += coulomb_harmonic_CM(nn, n, l, hw)
    t2 = time.time()
    profiler.add_timing("From LSJ to HO", t2 - t1)
    return temp  # MeV: both momentum integrations use MeV, so no hc factor.


# potential in lab frame under HO basis,
# T = -1 for pp, 0 for pn, +1 for nn.
def v_harmonic_lab(
    nap: int,
    lap: int,
    jap: float,
    nbp: int,
    lbp: int,
    jbp: float,
    na: int,
    la: int,
    ja: float,
    nb: int,
    lb: int,
    jb: float,
    J: int,
    T: int,
):
    t1 = time.time()
    # check parity conservation first.
    if (lap + lbp + la + lb) % 2:
        return 0
    tan_beta = 1.0
    sa, sb, sap, sbp = 0.5, 0.5, 0.5, 0.5  # all spin one-half particles.
    norm_apbp = 1 / np.sqrt(1 + bm.dirac_delta_int(nap, nbp) * bm.dirac_delta_int(lap, lbp) * bm.dirac_delta_float(jap, jbp) * np.abs(T))
    norm_ab = 1 / np.sqrt(1 + bm.dirac_delta_int(na, nb) * bm.dirac_delta_int(la, lb) * bm.dirac_delta_float(ja, jb) * np.abs(T))
    norm_spe = 1 + np.abs(T)
    fac1 = norm_apbp * norm_ab * norm_spe  # normalization factor
    fac2 = bm.hat(ja) * bm.hat(jb) * bm.hat(jap) * bm.hat(jbp)
    E_max = 2 * na + la + 2 * nb + lb
    Ep_max = 2 * nap + lap + 2 * nbp + lbp
    temp = 0.0
    temp3 = 0
    for S in [0, 1]:
        L_min = np.abs(la - lb)
        L_max = la + lb
        Lp_min = np.abs(lap - lbp)
        Lp_max = lap + lbp
        for L in range(L_min, L_max + 1, 1):
            ninej_factor1 = f9j(la, sa, ja, lb, sb, jb, L, S, J)
            if np.abs(ninej_factor1) < 1e-8:
                continue
            if not is_triangle(L, S, J):
                continue
            for Lp in range(Lp_min, Lp_max + 1, 1):
                ninej_factor2 = f9j(lap, sap, jap, lbp, sbp, jbp, Lp, S, J)
                if np.abs(ninej_factor2) < 1e-8:
                    continue
                if not is_triangle(Lp, S, J):
                    continue
                temp2 = 0
                for Lambda in range(0, min(E_max, Ep_max) + 1, 1):
                    twoN_max = min(np.abs(E_max - Lambda), np.abs(Ep_max - Lambda))
                    for twoN in range(0, twoN_max + 1, 2):
                        N = int(twoN / 2)
                        lambda_max = E_max - 2 * N - Lambda
                        if lambda_max < 0:
                            continue
                        for lambdaa in range(0, lambda_max + 1, 1):
                            twon = E_max - 2 * N - Lambda - lambdaa
                            if (twon % 2) != 0 or (twon < 0):
                                continue
                            n = int(twon / 2)
                            lambdap_max = Ep_max - 2 * N - Lambda
                            if lambdap_max < 0:
                                continue
                            for lambdap in range(0, lambdap_max + 1, 1):
                                twonp = Ep_max - 2 * N - Lambda - lambdap
                                if (twonp % 2) != 0 or (twonp < 0):
                                    continue
                                npp = int(twonp / 2)
                                temp1 = 0
                                j_min = max(abs(lambdaa - S), abs(lambdap - S))
                                j_max = min(lambdaa + S, lambdap + S)
                                if j_min > j_max:
                                    continue
                                morsh1 = moshinsky(N, Lambda, n, lambdaa, na, la, nb, lb, L, tan_beta)
                                morsh2 = moshinsky(
                                    N,
                                    Lambda,
                                    npp,
                                    lambdap,
                                    nap,
                                    lap,
                                    nbp,
                                    lbp,
                                    Lp,
                                    tan_beta,
                                )
                                if abs(morsh1) < 1e-8 or abs(morsh2) < 1e-8:
                                    continue
                                if not is_triangle(Lambda, lambdaa, L):
                                    continue
                                if not is_triangle(Lambda, lambdap, Lp):
                                    continue
                                for j in range(j_min, j_max + 1, 1):
                                    if not is_triangle(j, Lambda, J):
                                        continue
                                    if j > JMAX:
                                        continue
                                    sixj_factor1 = f6j(Lambda, lambdaa, L, S, J, j)
                                    sixj_factor2 = f6j(Lambda, lambdap, Lp, S, J, j)
                                    if abs(sixj_factor1) < 1e-8 or abs(sixj_factor2) < 1e-8:
                                        continue
                                    temp_harmonic = v_harmonic_CM(npp, n, lambdap, lambdaa, S, j, T, alpha)
                                    temp1 = temp1 + sixj_factor1 * sixj_factor2 * (2.0 * j + 1.0) * temp_harmonic
                                temp2 = temp2 + temp1 * morsh1 * morsh2
                temp3 = temp3 + temp2 * ninej_factor1 * ninej_factor2 * (2.0 * S + 1.0) * (2.0 * L + 1.0) * (2.0 * Lp + 1.0)
    temp = temp3 * fac1 * fac2
    t2 = time.time()
    profiler.add_timing("From C.M. to Lab.", t2 - t1)
    return temp


# single particle basis,
# ( index , n , l , 2j , tz , e )
def gen_single_particle_basis(emax: int):
    basis_temp = []
    index = 1
    for e_temp in range(0, emax + 1, 1):
        for l in range(0, e_temp + 1, 1):
            for two_n in range(0, e_temp - l + 1, 1):
                n = int(two_n / 2)
                if (2 * n + l) != e_temp:
                    continue
                if l == 0:
                    twoj = 1
                    basis_temp.append((index, n, l, twoj, -1, e_temp))
                    basis_temp.append((index + 1, n, l, twoj, 1, e_temp))
                    index = index + 2
                else:
                    twoj_minus = 2 * l - 1
                    twoj_plus = 2 * l + 1
                    basis_temp.append((index, n, l, twoj_minus, -1, e_temp))
                    basis_temp.append((index + 1, n, l, twoj_minus, 1, e_temp))
                    basis_temp.append((index + 2, n, l, twoj_plus, -1, e_temp))
                    basis_temp.append((index + 3, n, l, twoj_plus, 1, e_temp))
                    index = index + 4
    return basis_temp


def gen_two_particle_basis(emax: int, e2max: int):
    basis = gen_single_particle_basis(emax)
    spnum = len(basis)
    temp = []
    for a in range(0, spnum, 1):
        for b in range(a, spnum, 1):
            idxa, na, la, twoja, tza, ea = basis[a]
            idxb, nb, lb, twojb, tzb, eb = basis[b]
            e2ab = ea + eb
            if e2ab > e2max:  # e2max truncation
                continue
            tzab = int((tza + tzb) / 2)
            parityab = 1 - 2 * ((la + lb) % 2)
            jabmin = np.abs(int((twoja - twojb) / 2))
            jabmax = np.abs(int((twoja + twojb) / 2))
            if tza == 1 and tzb == -1:  # Keep pn order for every J.
                idxa, idxb = idxb, idxa
            for jab in range(jabmin, jabmax + 1, 1):
                if a == b and jab % 2 == 1:  # jab must be even if a=b.
                    continue
                temp.append((idxa, idxb, jab, parityab, tzab, e2ab))
    return temp


def gen_channels(emax: int, e2max: int):
    temp = []
    tp_basis = gen_two_particle_basis(emax, e2max)
    tpnum = len(tp_basis)
    for s in range(0, tpnum, 1):
        for t in range(s, tpnum, 1):
            idxa, idxb, jab, parityab, tzab, e2ab = tp_basis[s]
            idxc, idxd, jcd, paritycd, tzcd, e2cd = tp_basis[t]
            if jab == jcd and parityab == paritycd and tzab == tzcd:
                temp.append([idxa, idxb, idxc, idxd, jab, 0])
    return temp


def write_sp_orbit_to_file(filename, arrays, width=6):
    num_columns = 5
    num_rows = len(arrays)
    formatted_rows = []

    for row_idx in range(num_rows):
        row_data = []
        for col_idx in range(num_columns):
            value = arrays[row_idx][col_idx]
            formatted_value = str(value)
            row_data.append(formatted_value.ljust(width))
        formatted_rows.append(" ".join(row_data))

    with open(filename, "w") as file:
        file.write(str(num_rows) + "\n")
        file.write("\n".join(formatted_rows))


def write_2bme_to_file(filename, arrays, width=6):
    num_columns = 6
    num_rows = len(arrays)
    formatted_rows = []

    for row_idx in range(num_rows):
        row_data = []
        for col_idx in range(num_columns):
            value = arrays[row_idx][col_idx]
            if col_idx == 5:  # float value at position 6
                formatted_value = format(value, ".8f")
            else:  # integer values at positions 1-5
                formatted_value = str(value)
            row_data.append(formatted_value.ljust(width))
        formatted_rows.append(" ".join(row_data))

    with open(filename, "w") as file:
        file.write("\n".join(formatted_rows))


################################################################################################################


################################################################################################################

utility.section_message("Generating Two-Body Channels")

t3 = time.time()

basis = gen_single_particle_basis(emax)
channels = gen_channels(emax, e2max)

t4 = time.time()
profiler.add_timing("Generating Two-Body Channels", t4 - t3)

################################################################################################################


################################################################################################################

utility.section_message("Calculating TBMEs")

with tqdm(total=len(channels), colour="green", disable=(not verbose)) as pbar:
    for channel in channels:
        idxa, idxb, idxc, idxd, jab, v_ini = channel
        idxa, na, la, twoja, tza, ea = basis[idxa - 1]
        idxb, nb, lb, twojb, tzb, eb = basis[idxb - 1]
        idxc, nc, lc, twojc, tzc, ec = basis[idxc - 1]
        idxd, nd, ld, twojd, tzd, ed = basis[idxd - 1]
        ja = twoja / 2.0
        jb = twojb / 2.0
        jc = twojc / 2.0
        jd = twojd / 2.0
        tzab = int((tza + tzb) / 2)
        v = v_harmonic_lab(na, la, ja, nb, lb, jb, nc, lc, jc, nd, ld, jd, jab, tzab)
        channel[5] = v
        # print(
        #     "    "
        #     + str(channel[0])
        #     + "    "
        #     + str(channel[1])
        #     + "    "
        #     + str(channel[2])
        #     + "    "
        #     + str(channel[3])
        #     + "    "
        #     + str(channel[4])
        #     + "    "
        #     + "%4.8f" % channel[5]
        # )
        pbar.update(1)
pbar.close()

################################################################################################################

################################################################################################################

utility.section_message("Storing TBMEs")

t5 = time.time()

sp_orbit_file = "./result/sp-orbit.snt"
coulomb_suffix = "-coulomb" if include_coulomb else ""
TBME_file = f"./result/TBME-{potential_type}-hw{hw}-emax{emax}-e2max{e2max}{coulomb_suffix}.snt"

# write_sp_orbit_to_file(sp_orbit_file, basis)
write_2bme_to_file(TBME_file, channels)

t6 = time.time()
profiler.add_timing("Storing TBMEs", t6 - t5)

################################################################################################################

utility.section_message("Timings")

profiler.print_timings()

################################################################################################################

utility.footer_message()
