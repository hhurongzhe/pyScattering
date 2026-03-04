import constants as const
import numpy as np
import WignerSymbol as ws
import scipy
import sys
import os
import time
import utility
import profiler
import basic_math as bm
import scipy.special as special_function


class nn_studio:
    def __init__(self, jmin, jmax, tz, Np=75, mesh_type="gauleg_infinite"):
        self.jmin = jmin
        self.jmax = jmax
        self.tz = tz
        self.channels_uncoupled = []  # (l, s, j, t, tz)
        self.channels_coupled = []  # (j, t, tz)
        self.setup_pw_channels()

        self.Np = Np
        if mesh_type == "gauleg_infinite":
            self.pmesh, self.wmesh = self.gauss_legendre_inf_mesh()
        elif mesh_type == "gauleg_finite":
            self.pmesh, self.wmesh = self.gauss_legendre_line_mesh(1e-16, 1200)

        self.Tlabs = None

        # potential
        self.V = None

        # Tmatrices
        self.Tmtx_uncoupled = []
        self.Tmtx_coupled = []
        self.Tmtx_onshell_uncoupled = []
        self.Tmtx_onshell_coupled = []

        # phase shifts
        self.phase_shifts_uncoupled = []
        self.phase_shifts_coupled = []

        # theta in degree
        self.theta = np.linspace(1, 180, 180)

        # spin matrixs
        self.m11 = []
        self.m10 = []
        self.mpm = []
        self.m01 = []
        self.m00 = []
        self.mss = []

        # spin observables
        self.spin_obs = {}

    def setup_pw_channels(self):
        if self.jmin == 0:
            self.channels_uncoupled.append((0, 0, 0, 1, self.tz))  # 1S0
            self.channels_uncoupled.append((1, 1, 0, 1, self.tz))  # 3P0
        # uncoupled channels
        for j in range(self.jmin, self.jmax + 1):
            if j == 0:
                continue
            for s in [0, 1]:
                t = 1 - (j + s) % 2
                if self.tz != 0 and t == 0:
                    continue
                self.channels_uncoupled.append((j, s, j, t, self.tz))
        # coupled channels
        for j in range(self.jmin, self.jmax + 1):
            t = 1 - j % 2
            if self.tz != 0 and t == 0:
                continue
            if j == 0:
                continue
            self.channels_coupled.append((j, t, self.tz))

    def gauss_legendre_line_mesh(self, a, b, num=None):
        if num is None:
            x, w = np.polynomial.legendre.leggauss(self.Np)
        else:
            x, w = np.polynomial.legendre.leggauss(num)
        # Translate x values from the interval [-1, 1] to [a, b]
        t = 0.5 * (x + 1) * (b - a) + a
        u = w * 0.5 * (b - a)
        return t, u

    def gauss_legendre_inf_mesh(self):
        scale = 100.0
        x, w = np.polynomial.legendre.leggauss(self.Np)
        pi_over_4 = np.pi / 4.0
        t = scale * np.tan(pi_over_4 * (x + 1.0))
        u = scale * pi_over_4 / np.cos(pi_over_4 * (x + 1.0)) ** 2 * w
        return t, u

    def tzname(self):
        if self.tz == -1:
            return "pp"
        elif self.tz == 0:
            return "np"
        elif self.tz == +1:
            return "nn"
        else:
            exit("unknown isospin projection")

    def lab2rel(self, Tlab):
        if self.tz == -1:
            mu = const.Mp / 2
            ko2 = 0.5 * const.Mp * Tlab
        elif self.tz == 0:
            mu = const.Mp * const.Mn / (const.Mp + const.Mn)
            ko2 = const.Mp**2 * Tlab * (Tlab + 2 * const.Mn) / ((const.Mp + const.Mn) ** 2 + 2 * Tlab * const.Mp)
        elif self.tz == +1:
            mu = const.Mn / 2
            ko2 = 0.5 * const.Mn * Tlab
        else:
            exit("unknown isospin projection")
        if ko2 < 0:
            ko = np.complex(0, np.sqrt(np.abs(ko2)))
        else:
            ko = np.sqrt(ko2)
        return ko, mu

    def Vmtx(self, this_mesh, ll, l, s, j, t, tz):
        mtx = np.zeros((len(this_mesh), len(this_mesh)))
        for pidx, p in enumerate(this_mesh):
            for ppidx, pp in enumerate(this_mesh):
                mtx[ppidx][pidx] = self.V.potential(ll, l, pp, p, j, s, tz)
        return np.array(mtx)

    def VCmtx(self, this_mesh, ll, l, s, j, t, tz, ko):
        mtx = np.zeros((len(this_mesh), len(this_mesh)))
        for pidx, p in enumerate(this_mesh):
            for ppidx, pp in enumerate(this_mesh):
                mtx[ppidx][pidx] = self.V.potential_cutoff_coulomb(ll, l, pp, p, j, s, tz, ko)
        return np.array(mtx)

    def setup_Vmtx(self, is_coupled, this_channel, ko=None):
        if ko is None:
            this_mesh = self.pmesh
        else:
            this_mesh = np.hstack((self.pmesh, ko))
        if not is_coupled:
            (l, s, j, t, tz) = this_channel
            V = np.copy(self.Vmtx(this_mesh, l, l, s, j, t, tz))
            if tz == -1:
                V += np.copy(self.VCmtx(this_mesh, l, l, s, j, t, tz, ko))
            return V
        else:
            (j, t, tz) = this_channel
            Vmm = np.copy(self.Vmtx(this_mesh, j - 1, j - 1, 1, j, t, tz))
            Vmp = np.copy(self.Vmtx(this_mesh, j - 1, j + 1, 1, j, t, tz))
            Vpm = np.copy(self.Vmtx(this_mesh, j + 1, j - 1, 1, j, t, tz))
            Vpp = np.copy(self.Vmtx(this_mesh, j + 1, j + 1, 1, j, t, tz))
            if tz == -1:
                Vmm += np.copy(self.VCmtx(this_mesh, j - 1, j - 1, 1, j, t, tz, ko))
                Vpp += np.copy(self.VCmtx(this_mesh, j + 1, j + 1, 1, j, t, tz, ko))
            V = np.copy(np.vstack((np.hstack((Vmm, Vmp)), np.hstack((Vpm, Vpp)))))
            return V

    def setup_G0_vector(self, is_coupled, ko, mu):
        G0_vec = np.zeros((self.Np + 1), dtype=complex)
        G0_vec[0 : self.Np] = 2 * mu * self.wmesh * self.pmesh**2 / (ko**2 - self.pmesh**2)
        G0_vec[self.Np] = -2 * mu * np.sum(self.wmesh / (ko**2 - self.pmesh**2)) * ko**2 - 1j * np.pi * ko * mu
        return G0_vec

    def setup_VG0_kernel(self, is_coupled, Vmtx, ko, mu):
        G0_vec = self.setup_G0_vector(is_coupled, ko, mu)
        VG0 = np.zeros(Vmtx.shape, dtype=complex)
        if not is_coupled:
            for i in range(self.Np + 1):
                VG0[:, i] = Vmtx[:, i] * G0_vec[i]
        else:
            for i in range(2 * (self.Np + 1)):
                VG0[:, i] = Vmtx[:, i] * G0_vec[i % (self.Np + 1)]
        return VG0

    def solve_lippmann_schwinger(self, is_coupled, Vmtx, ko, mu):
        # matrix inversion: T = V + VGT --> T = (1-VG)^{-1}V
        VG0 = self.setup_VG0_kernel(is_coupled, Vmtx, ko, mu)
        if not is_coupled:
            dim = self.Np + 1
        else:
            dim = 2 * (self.Np + 1)
        Id = np.identity(dim)
        T = np.linalg.solve(Id - VG0, Vmtx)
        return T

    def compute_phase_shift_uncoupled(self, ko, mu, T):
        rad2deg = 180.0 / np.pi
        fac = np.pi * mu * ko
        S = 1 - 2j * fac * T
        if np.abs(np.abs(S) - 1.0) > 1e-6:
            sys.exit("Error: S-matrix is not unitary in uncoupled channel!")
        delta = (-0.5 * 1j) * np.log(S)
        return np.real(delta * rad2deg)

    def compute_phase_shift_coupled(self, ko, mu, T11, T12, T22):
        rad2deg = 180.0 / np.pi
        fac = np.pi * mu * ko
        T = np.array([[T11, T12], [T12, T22]])
        S = np.identity(2) - 2j * fac * T
        if np.abs(np.abs(np.linalg.det(S)) - 1.0) > 1e-6:
            sys.exit("Error: S-matrix is not unitary in coupled channel!")
        # Blatt-Biedenharn (BB) convention
        twoEpsilonJ_BB = np.arctan(2 * T12 / (T11 - T22))  # mixing parameter
        delta_plus_BB = -0.5 * 1j * np.log(1 - 1j * fac * (T11 + T22) + 1j * fac * (2 * T12) / np.sin(twoEpsilonJ_BB))
        delta_minus_BB = -0.5 * 1j * np.log(1 - 1j * fac * (T11 + T22) - 1j * fac * (2 * T12) / np.sin(twoEpsilonJ_BB))
        # Stapp convention (bar-phase shifts)
        cos2e = np.cos(twoEpsilonJ_BB / 2) * np.cos(twoEpsilonJ_BB / 2)
        cos_2dp = np.cos(2 * delta_plus_BB)
        cos_2dm = np.cos(2 * delta_minus_BB)
        sin_2dp = np.sin(2 * delta_plus_BB)
        sin_2dm = np.sin(2 * delta_minus_BB)
        aR = np.real(cos2e * cos_2dm + (1 - cos2e) * cos_2dp)
        aI = np.real(cos2e * sin_2dm + (1 - cos2e) * sin_2dp)
        delta_minus = 0.5 * np.arctan2(aI, aR)
        aR = np.real(cos2e * cos_2dp + (1 - cos2e) * cos_2dm)
        aI = np.real(cos2e * sin_2dp + (1 - cos2e) * sin_2dm)
        delta_plus = 0.5 * np.arctan2(aI, aR)
        tmp = 0.5 * np.sin(twoEpsilonJ_BB)
        aR = tmp * (cos_2dm - cos_2dp)
        aI = tmp * (sin_2dm - sin_2dp)
        tmp = delta_plus + delta_minus
        epsilon = 0.5 * np.arcsin(aI * np.cos(tmp) - aR * np.sin(tmp))
        if ko < 150:
            if delta_minus < 0:
                delta_minus += np.pi
                epsilon *= -1.0
        return [np.real(delta_minus * rad2deg), np.real(delta_plus * rad2deg), np.real(epsilon * rad2deg)]

    def compute_phase_shift_uncoupled_matched(self, ko, mu, Ts, this_channel):
        rad2deg = 180.0 / np.pi
        fac = np.pi * mu * ko
        Rc = 10.0 / const.hbarc
        z = ko * Rc
        factor_rel = (1 + 2 * ko * ko / const.Mp**2) / (np.sqrt(1 + ko * ko / const.Mp**2))
        eta = mu / ko * const.alpha * factor_rel
        (l, s, j, t, tz) = this_channel
        S = 1 - 2j * fac * Ts
        deltaS = (-0.5 * 1j) * np.log(S)
        # matching...
        AL0 = (bm.F(l, 0, z) + bm.G(l, 0, z) * np.tan(deltaS)) / (bm.dF(l, 0, z) + bm.dG(l, 0, z) * np.tan(deltaS))
        deltaC = np.arctan((AL0 * bm.dF(l, eta, z) - bm.F(l, eta, z)) / (bm.G(l, eta, z) - AL0 * bm.dG(l, eta, z)))
        return np.real(deltaC * rad2deg)

    def compute_phase_shift_coupled_matched(self, ko, mu, Ts11, Ts12, Ts22, this_channel):
        rad2deg = 180.0 / np.pi
        fac = np.pi * mu * ko
        Rc = 10.0 / const.hbarc
        z = ko * Rc
        factor_rel = (1 + 2 * ko * ko / const.Mp**2) / (np.sqrt(1 + ko * ko / const.Mp**2))
        eta = mu / ko * const.alpha * factor_rel
        # matching...
        Ts = -fac * np.array([[Ts11, Ts12], [Ts12, Ts22]])
        A = np.identity(2) + 1j * Ts
        Rs = Ts @ np.linalg.inv(A)
        (j, t, tz) = this_channel
        jm = j - 1
        jp = j + 1
        F0mtx = np.array([[bm.F(jm, 0, z), 0], [0, bm.F(jp, 0, z)]]).astype(np.float64)
        dF0mtx = np.array([[bm.dF(jm, 0, z), 0], [0, bm.dF(jp, 0, z)]]).astype(np.float64)
        G0mtx = np.array([[bm.G(jm, 0, z), 0], [0, bm.G(jp, 0, z)]]).astype(np.float64)
        dG0mtx = np.array([[bm.dG(jm, 0, z), 0], [0, bm.dG(jp, 0, z)]]).astype(np.float64)
        F1mtx = np.array([[bm.F(jm, eta, z), 0], [0, bm.F(jp, eta, z)]]).astype(np.float64)
        dF1mtx = np.array([[bm.dF(jm, eta, z), 0], [0, bm.dF(jp, eta, z)]]).astype(np.float64)
        G1mtx = np.array([[bm.G(jm, eta, z), 0], [0, bm.G(jp, eta, z)]]).astype(np.float64)
        dG1mtx = np.array([[bm.dG(jm, eta, z), 0], [0, bm.dG(jp, eta, z)]]).astype(np.float64)
        A0 = (F0mtx + G0mtx @ Rs) @ np.linalg.inv(dF0mtx + dG0mtx @ Rs)
        Rc = np.linalg.inv(G1mtx - A0 @ dG1mtx) @ (A0 @ dF1mtx - F1mtx)
        Tc = -(1 / fac) * np.linalg.inv(np.identity(2) - 1j * Rc) @ Rc
        T11c, T12c, T22c = Tc[0, 0], Tc[0, 1], Tc[1, 1]
        # get matched phase shifts
        twoEpsilonJ_BB = np.arctan(2 * T12c / (T11c - T22c))
        delta_plus_BB = -0.5 * 1j * np.log(1 - 1j * fac * (T11c + T22c) + 1j * fac * (2 * T12c) / np.sin(twoEpsilonJ_BB))
        delta_minus_BB = -0.5 * 1j * np.log(1 - 1j * fac * (T11c + T22c) - 1j * fac * (2 * T12c) / np.sin(twoEpsilonJ_BB))
        # Stapp convention (bar-phase shifts)
        cos2e = np.cos(twoEpsilonJ_BB / 2) * np.cos(twoEpsilonJ_BB / 2)
        cos_2dp = np.cos(2 * delta_plus_BB)
        cos_2dm = np.cos(2 * delta_minus_BB)
        sin_2dp = np.sin(2 * delta_plus_BB)
        sin_2dm = np.sin(2 * delta_minus_BB)
        aR = np.real(cos2e * cos_2dm + (1 - cos2e) * cos_2dp)
        aI = np.real(cos2e * sin_2dm + (1 - cos2e) * sin_2dp)
        delta_minus = 0.5 * np.arctan2(aI, aR)
        aR = np.real(cos2e * cos_2dp + (1 - cos2e) * cos_2dm)
        aI = np.real(cos2e * sin_2dp + (1 - cos2e) * sin_2dm)
        delta_plus = 0.5 * np.arctan2(aI, aR)
        tmp = 0.5 * np.sin(twoEpsilonJ_BB)
        aR = tmp * (cos_2dm - cos_2dp)
        aI = tmp * (sin_2dm - sin_2dp)
        tmp = delta_plus + delta_minus
        epsilon = 0.5 * np.arcsin(aI * np.cos(tmp) - aR * np.sin(tmp))
        return [np.real(delta_minus * rad2deg), np.real(delta_plus * rad2deg), np.real(epsilon * rad2deg)]

    def compute_Tmtx(self):
        for is_coupled, channels in [(False, self.channels_uncoupled), (True, self.channels_coupled)]:
            for channel in channels:
                if not is_coupled:
                    (l, s, j, t, tz) = channel
                else:
                    (j, t, tz) = channel
                phase_shifts_for_this_channel = []
                Tmtx_onshell_for_this_channel = []
                for Tlab in self.Tlabs:
                    print(f"Tlab = {Tlab} MeV")
                    ko, mu = self.lab2rel(Tlab)
                    t1 = time.time()
                    Vmtx = self.setup_Vmtx(is_coupled, channel, ko)
                    t2 = time.time()
                    profiler.add_timing("Setup V Matrix", t2 - t1)
                    this_T = self.solve_lippmann_schwinger(is_coupled, Vmtx, ko, mu)
                    t3 = time.time()
                    profiler.add_timing("Solve LS", t3 - t2)
                    if not is_coupled:
                        self.Tmtx_uncoupled.append(this_T)
                        T_on_shell = this_T[-1, -1]
                        if tz == -1:
                            this_phase_shift = self.compute_phase_shift_uncoupled_matched(ko, mu, T_on_shell, channel)
                        else:
                            this_phase_shift = self.compute_phase_shift_uncoupled(ko, mu, T_on_shell)
                    else:
                        self.Tmtx_coupled.append(this_T)
                        Np = int((this_T.shape[0] - 2) / 2)
                        Tmm = this_T[Np, Np]
                        Tmp = this_T[Np, 2 * Np + 1]
                        Tpp = this_T[2 * Np + 1, 2 * Np + 1]
                        T_on_shell = [Tmm, Tmp, Tpp]
                        if tz == -1:
                            this_phase_shift = self.compute_phase_shift_coupled_matched(ko, mu, Tmm, Tmp, Tpp, channel)
                        else:
                            this_phase_shift = self.compute_phase_shift_coupled(ko, mu, Tmm, Tmp, Tpp)
                    t4 = time.time()
                    profiler.add_timing("Solve Phase Shifts", t4 - t3)
                    Tmtx_onshell_for_this_channel.append(T_on_shell)
                    phase_shifts_for_this_channel.append(this_phase_shift)
                if not is_coupled:
                    self.Tmtx_onshell_uncoupled.append(np.array(Tmtx_onshell_for_this_channel))
                    self.phase_shifts_uncoupled.append(np.array(phase_shifts_for_this_channel))
                else:
                    self.Tmtx_onshell_coupled.append(np.array(Tmtx_onshell_for_this_channel))
                    self.phase_shifts_coupled.append(np.array(phase_shifts_for_this_channel))

    @staticmethod
    def sph_harm(l, m, theta, phi):
        return scipy.special.sph_harm_y(l, m, theta, phi)

    def get_spin_coeff(self, Sp, S, mp, m, the, lp, l, j):
        fac1 = np.sqrt(4 * np.pi * (2 * l + 1))
        fac2 = (1j) ** (l - lp)
        fac3 = ws.CG(2 * l, 2 * S, 2 * j, 0, 2 * m, 2 * m)
        fac4 = ws.CG(2 * lp, 2 * Sp, 2 * j, 2 * (m - mp), 2 * mp, 2 * m)
        fac5 = self.sph_harm(lp, m - mp, the, 0)
        fac = fac1 * fac2 * fac3 * fac4 * fac5
        return fac

    def get_mmatrix_coulomb(self, the, koid):
        ko, mu = self.lab2rel(self.Tlabs[koid])
        factor_rel = (1 + 2 * ko * ko / const.Mp**2) / np.sqrt(1 + ko * ko / const.Mp**2)
        eta = mu / ko * const.alpha * factor_rel
        temp = -1j * eta * np.log(np.sin(the / 2) ** 2)
        Mc = -eta / (2 * ko * np.sin(the / 2) ** 2) * np.exp(temp)
        return Mc

    @staticmethod
    def get_coulomb_phase(l, eta):
        return np.angle(special_function.gamma(l + 1 + 1j * eta))

    def compute_mmatrix(self, S, mp, m, the, koid):
        ko, mu = self.lab2rel(self.Tlabs[koid])
        # For pp, reconstruct nuclear amplitudes from Coulomb-matched phase shifts
        # and include Coulomb phase factors to keep the pp formula chain consistent.
        if self.tz == -1:
            factor_rel = (1 + 2 * ko * ko / const.Mp**2) / np.sqrt(1 + ko * ko / const.Mp**2)
            eta = mu / ko * const.alpha * factor_rel
            sigma_cache = {}

            def sigma(l):
                if l not in sigma_cache:
                    sigma_cache[l] = self.get_coulomb_phase(l, eta)
                return sigma_cache[l]

            temp = 0j
            for chan_idx, channel in enumerate(self.channels_uncoupled):
                (l, s, j, t, tz) = channel
                if s != S:
                    continue
                delta = np.deg2rad(self.phase_shifts_uncoupled[chan_idx][koid])
                amp = np.exp(2j * sigma(l)) * (np.exp(2j * delta) - 1) / (2j * ko)
                fac = self.get_spin_coeff(s, s, mp, m, the, l, l, j)
                temp += fac * amp

            if S == 1:
                for chan_idx, channel in enumerate(self.channels_coupled):
                    (j, t, tz) = channel
                    dm, dp, de = np.deg2rad(self.phase_shifts_coupled[chan_idx][koid])
                    sm = sigma(j - 1)
                    sp = sigma(j + 1)
                    amp_m = np.exp(2j * sm) * (np.cos(2 * de) * np.exp(2j * dm) - 1) / (2j * ko)
                    amp_p = np.exp(2j * sp) * (np.cos(2 * de) * np.exp(2j * dp) - 1) / (2j * ko)
                    amp_e = np.exp(1j * (sm + sp)) * (1j * np.sin(2 * de) * np.exp(1j * (dm + dp))) / (2j * ko)
                    temp += self.get_spin_coeff(1, 1, mp, m, the, j - 1, j - 1, j) * amp_m
                    temp += self.get_spin_coeff(1, 1, mp, m, the, j - 1, j + 1, j) * amp_e
                    temp += self.get_spin_coeff(1, 1, mp, m, the, j + 1, j - 1, j) * amp_e
                    temp += self.get_spin_coeff(1, 1, mp, m, the, j + 1, j + 1, j) * amp_p

            return 2 * temp

        temp = 0
        for chan_idx, channel in enumerate(self.channels_uncoupled):
            (l, s, j, t, tz) = channel
            if s != S:
                continue
            Tmtx_onshell_this_channel = self.Tmtx_onshell_uncoupled[chan_idx]
            Tmtx_onshell_this_ko = Tmtx_onshell_this_channel[koid]
            T = Tmtx_onshell_this_ko
            fac = self.get_spin_coeff(s, s, mp, m, the, l, l, j)
            temp += fac * T
        if S == 1:
            for chan_idx, channel in enumerate(self.channels_coupled):
                (j, t, tz) = channel
                Tmtx_onshell_this_channel = self.Tmtx_onshell_coupled[chan_idx]
                Tmtx_onshell_this_ko = Tmtx_onshell_this_channel[koid]
                Tmm, Tmp, Tpp = Tmtx_onshell_this_ko
                temp += self.get_spin_coeff(1, 1, mp, m, the, j - 1, j - 1, j) * Tmm
                temp += self.get_spin_coeff(1, 1, mp, m, the, j - 1, j + 1, j) * Tmp
                temp += self.get_spin_coeff(1, 1, mp, m, the, j + 1, j - 1, j) * Tmp
                temp += self.get_spin_coeff(1, 1, mp, m, the, j + 1, j + 1, j) * Tpp
        return -np.pi * mu * temp

    def get_all_mmatrix(self, the, koid):
        m11 = self.compute_mmatrix(1, 1, 1, the, koid)
        m10 = self.compute_mmatrix(1, 1, 0, the, koid)
        mpm = self.compute_mmatrix(1, 1, -1, the, koid)
        m01 = self.compute_mmatrix(1, 0, 1, the, koid)
        m00 = self.compute_mmatrix(1, 0, 0, the, koid)
        mss = self.compute_mmatrix(0, 0, 0, the, koid)
        check = m11 - m00 - np.sqrt(2) * (m10 + m01) * np.cos(the) / np.sin(the) - mpm
        if np.abs(check) > 1e-9:
            sys.exit(f"Error: Spin matrix relation not satisfied: {check}")
        return m11, m10, mpm, m01, m00, mss

    def build_m_matrix(self):
        ws.init(20, "Jmax", 3)
        self.m11 = []
        self.m10 = []
        self.mpm = []
        self.m01 = []
        self.m00 = []
        self.mss = []
        for idx, Tlab in enumerate(self.Tlabs):
            temp_m11 = []
            temp_m10 = []
            temp_mpm = []
            temp_m01 = []
            temp_m00 = []
            temp_mss = []
            for the in self.theta:
                # from degree to rad
                the = the * np.pi / 180
                m11, m10, mpm, m01, m00, mss = self.get_all_mmatrix(the, idx)
                temp_m11.append(m11)
                temp_m10.append(m10)
                temp_mpm.append(mpm)
                temp_m01.append(m01)
                temp_m00.append(m00)
                temp_mss.append(mss)
            self.m11.append(temp_m11)
            self.m10.append(temp_m10)
            self.mpm.append(temp_mpm)
            self.m01.append(temp_m01)
            self.m00.append(temp_m00)
            self.mss.append(temp_mss)

    def cal_spin_observables(self):
        self.spin_obs = {}
        self.spin_obs["DSG"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["D"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["P"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["A"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["R"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["Rp"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["Axx"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["Azz"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["Axz"] = np.zeros((len(self.Tlabs), len(self.theta)))
        self.spin_obs["TSG"] = np.zeros(len(self.Tlabs))
        for idx, Tlab in enumerate(self.Tlabs):
            for idx_the, the in enumerate(self.theta):
                the = the * np.pi / 180
                m11 = self.m11[idx][idx_the]
                m10 = self.m10[idx][idx_the]
                mpm = self.mpm[idx][idx_the]
                m01 = self.m01[idx][idx_the]
                m00 = self.m00[idx][idx_the]
                mss = self.mss[idx][idx_the]
                if self.tz == -1:
                    fcs = self.get_mmatrix_coulomb(the, idx) + self.get_mmatrix_coulomb(np.pi - the, idx)
                    fca = self.get_mmatrix_coulomb(the, idx) - self.get_mmatrix_coulomb(np.pi - the, idx)
                    mss += fcs
                    m11 += fca
                    m00 += fca
                I0 = 0.5 * np.abs(m11) ** 2 + 0.5 * np.abs(m10) ** 2 + 0.5 * np.abs(mpm) ** 2 + 0.5 * np.abs(m01) ** 2 + 0.25 * np.abs(m00) ** 2 + 0.25 * np.abs(mss) ** 2
                I01D = 0.25 * np.abs(m11 + mpm - mss) ** 2 + 0.25 * np.abs(m11 - mpm - m00) ** 2 + 0.5 * np.abs(m10 + m01) ** 2
                I0P = np.sqrt(2) / 4.0 * np.real(1j * (m10 - m01) * np.conj(m11 - mpm + m00))
                I0A = -0.5 * np.real((m00 + np.sqrt(2) * (np.cos(the) + 1) / np.sin(the) * m10) * np.conj(m11 + mpm + mss) - np.sqrt(2) / np.sin(the) * (m10 + m01) * np.conj(m11 + mpm)) * np.sin(the / 2.0)
                I0R = 0.5 * np.real((m00 + np.sqrt(2) * (np.cos(the) - 1) / np.sin(the) * m10) * np.conj(m11 + mpm + mss) + np.sqrt(2) / np.sin(the) * (m10 + m01) * np.conj(mss)) * np.cos(the / 2)
                I0Rp = 0.5 * np.real((m00 + np.sqrt(2) * (np.cos(the) + 1) / np.sin(the) * m10) * np.conj(m11 + mpm + mss) - np.sqrt(2) / np.sin(the) * (m10 + m01) * np.conj(mss)) * np.sin(the / 2)
                I0Axx = 0.25 * np.abs(m00) ** 2 - 0.25 * np.abs(mss) ** 2 - 0.5 * np.abs(m01) ** 2 + 0.5 * np.abs(m10) ** 2 + np.real(m11 * np.conj(mpm))
                I0Azz = 0.5 * np.abs(m11) ** 2 - 0.25 * np.abs(m00) ** 2 - 0.25 * np.abs(mss) ** 2 + 0.5 * np.abs(m01) ** 2 - 0.5 * np.abs(m10) ** 2 + 0.5 * np.abs(mpm) ** 2
                I0Axz = 0.25 * np.tan(the) * (np.abs(m11 - mpm) ** 2 - np.abs(m00) ** 2) - 0.5 / np.tan(the) * (np.abs(m01) ** 2 - np.abs(m10) ** 2)
                self.spin_obs["DSG"][idx][idx_the] = I0 * 10 * const.hbarc**2
                self.spin_obs["D"][idx][idx_the] = 1 - I01D / I0
                self.spin_obs["P"][idx][idx_the] = I0P / I0
                self.spin_obs["A"][idx][idx_the] = I0A / I0
                self.spin_obs["R"][idx][idx_the] = I0R / I0
                self.spin_obs["Rp"][idx][idx_the] = I0Rp / I0
                self.spin_obs["Axx"][idx][idx_the] = I0Axx / I0
                self.spin_obs["Azz"][idx][idx_the] = I0Azz / I0
                self.spin_obs["Axz"][idx][idx_the] = I0Axz / I0
            TSG_this = np.sum(self.spin_obs["DSG"][idx] * np.sin(self.theta * np.pi / 180)) * np.pi / 180
            self.spin_obs["TSG"][idx] = TSG_this

    @staticmethod
    def write_arrays_to_file(filename, column_names, arrays, width=24, precision=4):
        num_columns = len(arrays)
        num_rows = len(arrays[0])
        formatted_rows = []

        for row_idx in range(num_rows):
            row_data = []
            for col_idx in range(num_columns):
                value = arrays[col_idx][row_idx]
                formatted_value = format(value, f".{precision}f")
                row_data.append(formatted_value.ljust(width))
            formatted_rows.append(" ".join(row_data))

        with open(filename, "w") as file:
            names_line = " ".join(name.ljust(width) for name in column_names) + "\n"
            file.write(names_line)
            file.write("\n".join(formatted_rows))

    # Writter of observables
    def store_observables(self):
        dir = "result"  # observables files are generated in this dir.
        if not os.path.exists(dir):
            os.makedirs(dir)
        for idx, tlab in enumerate(self.Tlabs):
            file_name = "./" + dir + "/" + f"spin_observables_{self.V.chiral_type}_tlab{tlab}_Jmax{self.jmax}_{self.tzname()}.txt"
            name_list = ["theta", "DSG", "D", "P", "A", "R", "Rp", "Axx", "Azz", "Axz"]
            self.write_arrays_to_file(
                file_name,
                name_list,
                [
                    self.theta,
                    self.spin_obs["DSG"][idx],
                    self.spin_obs["D"][idx],
                    self.spin_obs["P"][idx],
                    self.spin_obs["A"][idx],
                    self.spin_obs["R"][idx],
                    self.spin_obs["Rp"][idx],
                    self.spin_obs["Axx"][idx],
                    self.spin_obs["Azz"][idx],
                    self.spin_obs["Axz"][idx],
                ],
            )

    def make_partial_wave_name(self, is_coupled, channel):
        orbit_table = "SPDFGHIJKLMNOQRTUVWXYZ"
        if not is_coupled:
            (l, s, j, t, tz) = channel
            name = [f"{2*s+1}{orbit_table[l]}{j}"]
        else:
            (j, t, tz) = channel
            name_mm = f"3{orbit_table[j-1]}{j}"
            name_pp = f"3{orbit_table[j+1]}{j}"
            name_e = f"E{j}"
            name = [name_mm, name_pp, name_e]
        return name

    # Writter of phase shifts
    def store_phase_shifts(self):
        dir = "result"  # observables files are generated in this dir.
        if not os.path.exists(dir):
            os.makedirs(dir)
        file_name = "./" + dir + "/" + f"phase_shifts_{self.V.chiral_type}_Jmax{self.jmax}_{self.tzname()}.txt"
        name_list = ["Tlab"]
        phases = [self.Tlabs]
        # 1. Process Uncoupled Channels
        for idx, channel in enumerate(self.channels_uncoupled):
            name_list.extend(self.make_partial_wave_name(False, channel))
            phases.append(self.phase_shifts_uncoupled[idx])
        # 2. Process Coupled Channels
        for idx, channel in enumerate(self.channels_coupled):
            name_list.extend(self.make_partial_wave_name(True, channel))
            coupled_data = self.phase_shifts_coupled[idx]
            phases.append(coupled_data[:, 0])
            phases.append(coupled_data[:, 1])
            phases.append(coupled_data[:, 2])
        # 3. Write to file
        self.write_arrays_to_file(file_name, name_list, phases)
