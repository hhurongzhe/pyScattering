import sys

sys.path.append("./deps")
import os
import numpy as np
import struct
import utility
import profiler
import time
from scipy.special import spherical_jn
from scipy.special import gamma
import scipy.special as special_function
import constants as const
import force_module as pot
from functools import lru_cache


# two-nucleon potentials are defined here, most of them are chiral potentials, including:
# loemn450, loemn500, loemn550, nloemn450, nloemn500, nloemn550,
# n2loemn450, n2loemn500, n2loemn550, n3loemn450, n3loemn500, n3loemn550,
# n4loemn450, n4loemn500, n4loemn550,
# n3loem,
# idaholocal1, idaholocal2, idaholocal3, ..., idaholocal12,
# n2loopt, n2losat, cdbonn, av18,
# scslo1fm, scslo1p1fm, scslo1p2fm, ..., scsn4lo1fm, scsn4lo1p1fm, scsn4lo1p2fm,
# losms400, losms450, losms500, losms550, ..., n4lo+sms400, n4lo+sms450, n4lo+sms500, n4lo+sms550.


# Chiral two-nucleon potential.
class two_nucleon_potential:
    def __init__(self, chiral_type, add_cutoff_coulomb=False):
        self.chiral_type = chiral_type
        # p mesh points and weights, in [MeV]
        self.pmesh_points, self.pmesh_weights = self.set_up_pmesh()
        # r mesh points and weights, in [fm]
        self.rmesh_points, self.rmesh_weights = self.set_up_rmesh()
        # if add cut-off coulomb potential, only useful when considering p-p scattering
        self.add_cutoff_coulomb = add_cutoff_coulomb

        # quadrature points for cos(theta) integration
        self.ntheta = 32
        # max j (affects precomputation of meshes for angular integrals)
        self.jmax = 10
        # z = cos(theta)
        self.z, self.w = np.polynomial.legendre.leggauss(self.ntheta)
        self.zJ = [self.z**J for J in range(0, self.jmax + 1)]
        self.P = []
        self.setup_legendre_polynomials()

    # potential in momentum space, p in [MeV] and V in [MeV^(-2)]
    @lru_cache(maxsize=None)
    def potential(self, ll, l, pp, p, j, s, tz):
        V = 0
        local_type = ["av18", "v8'", "idaholocal1", "idaholocal2", "idaholocal3", "idaholocal4", "idaholocal5", "idaholocal6", "idaholocal7", "idaholocal8", "idaholocal9", "idaholocal10", "idaholocal11", "idaholocal12"]
        if self.chiral_type in local_type:
            V = self.nonlocal_projection(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "loemn450":
            V = pot.loemn450(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "loemn500":
            V = pot.loemn500(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "loemn550":
            V = pot.loemn550(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "nloemn450":
            V = pot.nloemn450(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "nloemn500":
            V = pot.nloemn500(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "nloemn550":
            V = pot.nloemn550(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2loemn450":
            V = pot.n2loemn450(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2loemn500":
            V = pot.n2loemn500(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2loemn550":
            V = pot.n2loemn550(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n3loemn450":
            V = pot.n3loemn450(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n3loemn500":
            V = pot.n3loemn500(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n3loemn550":
            V = pot.n3loemn550(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4loemn450":
            V = pot.n4loemn450(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4loemn500":
            V = pot.n4loemn500(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4loemn550":
            V = pot.n4loemn550(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2loopt":
            V = pot.nnloopt(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2losat":
            V = pot.nnlosat(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "cdbonn":
            V = pot.cdbonnpot(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n3loem":
            V = pot.n3loem(ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scslo0p8fm":
            V = pot.scs(1, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scslo0p9fm":
            V = pot.scs(2, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scslo1fm":
            V = pot.scs(3, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scslo1p1fm":
            V = pot.scs(4, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scslo1p2fm":
            V = pot.scs(5, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsnlo0p8fm":
            V = pot.scs(1, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsnlo0p9fm":
            V = pot.scs(2, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsnlo1fm":
            V = pot.scs(3, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsnlo1p1fm":
            V = pot.scs(4, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsnlo1p2fm":
            V = pot.scs(5, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn2lo0p8fm":
            V = pot.scs(1, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn2lo0p9fm":
            V = pot.scs(2, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn2lo1fm":
            V = pot.scs(3, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn2lo1p1fm":
            V = pot.scs(4, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn2lo1p2fm":
            V = pot.scs(5, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn3lo0p8fm":
            V = pot.scs(1, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn3lo0p9fm":
            V = pot.scs(2, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn3lo1fm":
            V = pot.scs(3, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn3lo1p1fm":
            V = pot.scs(4, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn3lo1p2fm":
            V = pot.scs(5, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn4lo0p8fm":
            V = pot.scs(1, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn4lo0p9fm":
            V = pot.scs(2, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn4lo1fm":
            V = pot.scs(3, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn4lo1p1fm":
            V = pot.scs(4, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "scsn4lo1p2fm":
            V = pot.scs(5, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "losms400":
            V = pot.sms(1, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "losms450":
            V = pot.sms(2, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "losms500":
            V = pot.sms(3, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "losms550":
            V = pot.sms(4, 0, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "nlosms400":
            V = pot.sms(1, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "nlosms450":
            V = pot.sms(2, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "nlosms500":
            V = pot.sms(3, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "nlosms550":
            V = pot.sms(4, 1, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2losms400":
            V = pot.sms(1, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2losms450":
            V = pot.sms(2, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2losms500":
            V = pot.sms(3, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n2losms550":
            V = pot.sms(4, 2, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n3losms400":
            V = pot.sms(1, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n3losms450":
            V = pot.sms(2, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n3losms500":
            V = pot.sms(3, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n3losms550":
            V = pot.sms(4, 3, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4losms400":
            V = pot.sms(1, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4losms450":
            V = pot.sms(2, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4losms500":
            V = pot.sms(3, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4losms550":
            V = pot.sms(4, 4, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4lo+sms400":
            V = pot.sms(1, 5, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4lo+sms450":
            V = pot.sms(2, 5, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4lo+sms500":
            V = pot.sms(3, 5, ll, l, pp, p, j, s, tz)
        elif self.chiral_type == "n4lo+sms550":
            V = pot.sms(4, 5, ll, l, pp, p, j, s, tz)
        else:
            sys.exit(f"momentum space potential type {self.chiral_type} not implemented yet.")
        if self.add_cutoff_coulomb:
            V = V + self.potential_cutoff_coulomb(ll, l, pp, p, j, s, tz)
        return V

    # potential in position space, where r in [fm] and V in [MeV]
    @lru_cache(maxsize=None)
    def potential_local(self, ll, l, s, j, tz, r):
        V = 0
        nonlocal_type = [
            "loemn450",
            "loemn500",
            "loemn550",
            "nloemn450",
            "nloemn500",
            "nloemn550",
            "n2loemn450",
            "n2loemn500",
            "n2loemn550",
            "n3loemn450",
            "n3loemn500",
            "n3loemn550",
            "n4loemn450",
            "n4loemn500",
            "n4loemn550",
            "n2loopt",
            "n2losat",
            "cdbonn",
            "n3loem",
            "scslo0p8fm",
            "scslo0p9fm",
            "scslo1fm",
            "scslo1p1fm",
            "scslo1p2fm",
            "scsnlo0p8fm",
            "scsnlo0p9fm",
            "scsnlo1fm",
            "scsnlo1p1fm",
            "scsnlo1p2fm",
            "scsn2lo0p8fm",
            "scsn2lo0p9fm",
            "scsn2lo1fm",
            "scsn2lo1p1fm",
            "scsn2lo1p2fm",
            "scsn3lo0p8fm",
            "scsn3lo0p9fm",
            "scsn3lo1fm",
            "scsn3lo1p1fm",
            "scsn3lo1p2fm",
            "scsn4lo0p8fm",
            "scsn4lo0p9fm",
            "scsn4lo1fm",
            "scsn4lo1p1fm",
            "scsn4lo1p2fm",
            "losms400",
            "losms450",
            "losms500",
            "losms550",
            "nlosms400",
            "nlosms450",
            "nlosms500",
            "nlosms550",
            "n2losms400",
            "n2losms450",
            "n2losms500",
            "n2losms550",
            "n3losms400",
            "n3losms450",
            "n3losms500",
            "n3losms550",
            "n4losms400",
            "n4losms450",
            "n4losms500",
            "n4losms550",
            "n4lo+sms400",
            "n4lo+sms450",
            "n4lo+sms500",
            "n4lo+sms550",
        ]
        if self.chiral_type in nonlocal_type:
            V = self.local_projection(ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal1":
            type = 1
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal2":
            type = 2
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal3":
            type = 3
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal4":
            type = 4
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal5":
            type = 5
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal6":
            type = 6
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal7":
            type = 7
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal8":
            type = 8
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal9":
            type = 9
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal10":
            type = 10
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal11":
            type = 11
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "idaholocal12":
            type = 12
            V = pot.idahopot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "av18":
            type = 1
            V = pot.av18pot(type, ll, l, s, j, tz, r)
        elif self.chiral_type == "v8'":
            type = 2
            V = pot.av18pot(type, ll, l, s, j, tz, r)
        else:
            sys.exit(f"position space potential type {self.chiral_type} not implemented yet.")
        return V

    # Local projection based on Eq. (9) and (10) in K. A. Wendt, R. J. Furnstahl, and S. Ramanan, Phys. Rev. C 86, 014003 (2012).
    @lru_cache(maxsize=None)
    def local_projection(self, ll, l, s, j, tz, r):
        temp = 0
        if l == 0 or ll == 0:
            lll = np.max([l, ll])
            for idx, k in enumerate(self.pmesh_points):
                w = self.pmesh_weights[idx]
                temp = temp + (k**2) * w * spherical_jn(lll, k * r / const.hbarc) * self.potential(lll, 0, k, 0, j, s, tz)
            return temp
        else:
            norm = 4 / np.sqrt(np.pi) * gamma((l + 3) / 2) / gamma(l / 2)
            if ll != l:
                norm *= -1
            for idxkk, kk in enumerate(self.pmesh_points):
                wkk = self.pmesh_weights[idxkk]
                for idxk, k in enumerate(self.pmesh_points):
                    wk = self.pmesh_weights[idxk]
                    temp = temp + wkk * wk * (kk**2 / k) * spherical_jn(ll, kk * r / const.hbarc) * self.potential(ll, l, kk, k, j, s, tz)
            return temp * norm

    # Nonlocal projection based on Eq. (12) in K. A. Wendt, R. J. Furnstahl, and S. Ramanan, Phys. Rev. C 86, 014003 (2012).
    @lru_cache(maxsize=None)
    def nonlocal_projection(self, ll, l, pp, p, j, s, tz):
        temp = 0
        norm = 2 / np.pi
        for idx, r in enumerate(self.rmesh_points):
            w = self.rmesh_weights[idx]
            temp = temp + (r**2) * w * spherical_jn(l, p * r / const.hbarc) * spherical_jn(ll, pp * r / const.hbarc) * self.potential_local(ll, l, s, j, tz, r)
        return temp * norm / const.hbarc**3

    @staticmethod
    def legP(m, n, x):
        # [0]: get function value (not derivative)
        # [-1,-1]: get 'm', 'n' values
        return special_function.lpmn(m, n, x)[0][-1, -1]

    def setup_legendre_polynomials(self):
        fP = np.vectorize(self.legP, excluded={0, 1}, otypes=[np.float64])

        for j in range(0, self.jmax):
            this_P = fP(0, j, self.z)
            self.P.append(this_P)

    def pwd_integral(self, W, l, j):
        if j < 0:
            return
        return np.pi * np.sum(W * self.zJ[l] * self.P[j] * self.w)

    # partial-wave decomposition of central force
    def pwd_C(self, W, ll, l, pp, p, j, s):
        # uncoupled singlet
        if ll == l and l == j and s == 0:
            V = 2 * self.pwd_integral(W, 0, j)
            return V
        # uncoupled triplet
        elif ll == l and l == j and s == 1:
            V = 2 * self.pwd_integral(W, 0, j)
            return V
        elif ll != j and l != j and s == 1:
            if ll == (j + 1) and l == (j + 1):
                V = 2 * self.pwd_integral(W, 0, j + 1)
                return V
            # 3P0 case
            if j == 0:
                V = 2 * self.pwd_integral(W, 0, j + 1)
                return V
            if ll == (j - 1) and l == (j - 1):
                V = 2 * self.pwd_integral(W, 0, j - 1)
                return V
        return 0

    def potential_cutoff_coulomb(self, ll, l, pp, p, j, s, tz, ko):
        if tz != -1:
            return 0
        q2 = pp**2 + p**2 - 2 * pp * p * self.z
        epsilon = 1
        alpha = 1.0 / 137.035989
        factor = (1 + 2 * ko * ko / const.Mp**2) / np.sqrt(1 + ko * ko / const.Mp**2)
        factor_rela = np.sqrt(const.Mp / np.sqrt(const.Mp**2 + pp**2)) * np.sqrt(const.Mp / np.sqrt(const.Mp**2 + p**2))
        Cc = 4 * np.pi * epsilon * alpha * factor * factor_rela
        R = 10.0 / const.hbarc  # R = 10 fm
        Vc = Cc * (1 - np.cos(np.sqrt(q2) * R)) / q2
        V = self.pwd_C(Vc, ll, l, pp, p, j, s)
        # 1/(2pi)^3 normalization
        return (0.125 / np.pi**3) * V

    def gauss_legendre_line_mesh(self, a, b, Np):
        x, w = np.polynomial.legendre.leggauss(Np)
        # Translate x values from the interval [-1, 1] to [a, b]
        t = 0.5 * (x + 1) * (b - a) + a
        u = w * 0.5 * (b - a)
        return t, u

    def set_up_pmesh(self):
        n0 = 8
        mesh1p, mesh1pw = self.gauss_legendre_line_mesh(0.0, 200.0, 16 * n0)
        mesh2p, mesh2pw = self.gauss_legendre_line_mesh(200.0, 400.0, 8 * n0)
        mesh3p, mesh3pw = self.gauss_legendre_line_mesh(400.0, 800.0, 4 * n0)
        mesh4p, mesh4pw = self.gauss_legendre_line_mesh(800.0, 1600.0, 2 * n0)
        mesh5p, mesh5pw = self.gauss_legendre_line_mesh(1600.0, 3200.0, n0)
        pmesh_points = np.concatenate((mesh1p, mesh2p, mesh3p, mesh4p, mesh5p), axis=0)
        pmesh_weights = np.concatenate((mesh1pw, mesh2pw, mesh3pw, mesh4pw, mesh5pw), axis=0)
        return pmesh_points, pmesh_weights

    def set_up_rmesh(self):
        n0 = 4
        mesh1r, mesh1rw = self.gauss_legendre_line_mesh(0.0, 1.0, 64 * n0)
        mesh2r, mesh2rw = self.gauss_legendre_line_mesh(1.0, 2.0, 32 * n0)
        mesh3r, mesh3rw = self.gauss_legendre_line_mesh(2.0, 4.0, 16 * n0)
        mesh4r, mesh4rw = self.gauss_legendre_line_mesh(4.0, 8.0, 8 * n0)
        mesh5r, mesh5rw = self.gauss_legendre_line_mesh(8.0, 16.0, 4 * n0)
        mesh6r, mesh6rw = self.gauss_legendre_line_mesh(16.0, 32.0, 2 * n0)
        mesh7r, mesh7rw = self.gauss_legendre_line_mesh(32.0, 40.0, n0)
        rmesh_points = np.concatenate((mesh1r, mesh2r, mesh3r, mesh4r, mesh5r, mesh6r, mesh7r), axis=0)
        rmesh_weights = np.concatenate((mesh1rw, mesh2rw, mesh3rw, mesh4rw, mesh5rw, mesh6rw, mesh7rw), axis=0)
        # rmesh_points, rmesh_weights = self.gauss_legendre_line_mesh(0, 60, 600)
        return rmesh_points, rmesh_weights

    @staticmethod
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

    # Writter of interaction matrix elements in momentum space, where [kmax] = fm^(-1).
    # The resulting interaction mtx are in [MeV fm^(-3)].
    def write_mtx_momentum(self, kmax=8.0, N=100, Jmax=8):
        t1 = time.time()
        # this value of hw is in line with Miyagi's code,
        # used for unit transformation.
        hc = 197.32705
        hccubic = hc**3
        dir = "input_nn_files"  # files are generated in this dir.
        if not os.path.exists(dir):
            os.makedirs(dir)
        file_name = "./" + dir + "/" + f"{self.chiral_type}_kmax{kmax}_N{N}_Jmax{Jmax}.dat"
        channels = self.gen_mtx_channels(N, Jmax)
        MeshPoints, MeshWeights = self.gauss_legendre_line_mesh(0, kmax, N)
        with open(file_name, "w") as fp:
            fp.write(f"NMesh:\n{N}\n")
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
                if Ndim != N:
                    coup = True
                fp.write(f"J:\n{J}\n")
                fp.write(f"Prty:\n{Prty}\n")
                fp.write(f"S:\n{S}\n")
                fp.write(f"Tz:\n{Tz}\n")
                fp.write(f"Ndim:\n{Ndim}\n")
                fp.write("V:\n")
                if J == 0:
                    V = np.array([[self.potential(S, S, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints]) * hccubic
                    np.savetxt(fp, V, fmt="%.17f")
                else:
                    if not coup:
                        V = np.array([[self.potential(J, J, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints]) * hccubic
                        np.savetxt(fp, V, fmt="%.17f")
                    else:
                        Vpp = np.array([[self.potential(J + 1, J + 1, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints])
                        Vpm = np.array([[self.potential(J + 1, J - 1, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints])
                        Vmp = np.array([[self.potential(J - 1, J + 1, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints])
                        Vmm = np.array([[self.potential(J - 1, J - 1, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints])
                        V = np.block([[Vmm, Vmp], [Vpm, Vpp]]) * hccubic
                        np.savetxt(fp, V, fmt="%.17f")
        binary_file = "./" + dir + "/" + f"{self.chiral_type}_kmax{kmax}_N{N}_Jmax{Jmax}.bin"
        t2 = time.time()
        profiler.add_timing("Writing .dat", t2 - t1)
        self.convert_dat_to_binary(file_name, binary_file)
        t3 = time.time()
        profiler.add_timing("Writing .bin", t3 - t2)
        print("files are generated to :\n", file_name, "\n", binary_file)

    # same as write_mtx_momentum, but no CSB, no CIB
    def write_mtx_momentum_iso(self, kmax=8.0, N=100, Jmax=8):
        t1 = time.time()
        # this value of hw is in line with Miyagi's code,
        # used for unit transformation.
        hc = 197.32705
        hccubic = hc**3
        dir = "input_nn_files"  # files are generated in this dir.
        if not os.path.exists(dir):
            os.makedirs(dir)
        file_name = "./" + dir + "/" + f"{self.chiral_type}_iso_kmax{kmax}_N{N}_Jmax{Jmax}.dat"
        channels = self.gen_mtx_channels(N, Jmax)
        MeshPoints, MeshWeights = self.gauss_legendre_line_mesh(0, kmax, N)
        with open(file_name, "w") as fp:
            fp.write(f"NMesh:\n{N}\n")
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
                if Ndim != N:
                    coup = True
                fp.write(f"J:\n{J}\n")
                fp.write(f"Prty:\n{Prty}\n")
                fp.write(f"S:\n{S}\n")
                fp.write(f"Tz:\n{Tz}\n")
                fp.write(f"Ndim:\n{Ndim}\n")
                fp.write("V:\n")
                Tz = 0  # force Tz=0 for isospin symmetric potential
                if J == 0:
                    V = np.array([[self.potential(S, S, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints]) * hccubic
                    np.savetxt(fp, V, fmt="%.17f")
                else:
                    if not coup:
                        V = np.array([[self.potential(J, J, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints]) * hccubic
                        np.savetxt(fp, V, fmt="%.17f")
                    else:
                        Vpp = np.array([[self.potential(J + 1, J + 1, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints])
                        Vpm = np.array([[self.potential(J + 1, J - 1, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints])
                        Vmp = np.array([[self.potential(J - 1, J + 1, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints])
                        Vmm = np.array([[self.potential(J - 1, J - 1, pi, pj, J, S, Tz) for pj in MeshPoints] for pi in MeshPoints])
                        V = np.block([[Vmm, Vmp], [Vpm, Vpp]]) * hccubic
                        np.savetxt(fp, V, fmt="%.17f")
        binary_file = "./" + dir + "/" + f"{self.chiral_type}_iso_kmax{kmax}_N{N}_Jmax{Jmax}.bin"
        t2 = time.time()
        profiler.add_timing("Writing .dat", t2 - t1)
        self.convert_dat_to_binary(file_name, binary_file)
        t3 = time.time()
        profiler.add_timing("Writing .bin", t3 - t2)
        print("files are generated to :\n", file_name, "\n", binary_file)

    # Writter of interaction matrix elements in coordinate space, where [rmax] = fm.
    # The resulting interaction mtx are in [MeV].
    def write_mtx_coordinate(self, rmax=4.0, N=100, Jmax=8):
        t1 = time.time()
        # this value of hw is in line with Miyagi's code,
        # used for unit transformation.
        hc = 197.32705
        dir = "input_nn_files"  # files are generated in this dir.
        if not os.path.exists(dir):
            os.makedirs(dir)
        file_name = "./" + dir + "/" + f"{self.chiral_type}_rmax{rmax}_N{N}_Jmax{Jmax}.dat"
        channels = self.gen_mtx_channels(N, Jmax)
        MeshPoints, MeshWeights = self.gauss_legendre_line_mesh(0, rmax, N)
        with open(file_name, "w") as fp:
            fp.write(f"NMesh:\n{N}\n")
            fp.write(f"Jmax:\n{Jmax}\n")
            fp.write(f"NChan:\n{len(channels)}\n")
            fp.write("Coordinate Mesh Points:\n")
            np.savetxt(fp, MeshPoints, fmt="%.17f")
            fp.write("Coordinate Mesh Weights:\n")
            np.savetxt(fp, MeshWeights, fmt="%.17f")
            for chan in channels:
                J = chan[0]
                Prty = chan[1]
                S = chan[2]
                Tz = chan[3]
                Ndim = chan[4]
                coup = False
                if Ndim != N:
                    coup = True
                fp.write(f"J:\n{J}\n")
                fp.write(f"Prty:\n{Prty}\n")
                fp.write(f"S:\n{S}\n")
                fp.write(f"Tz:\n{Tz}\n")
                fp.write(f"Ndim:\n{Ndim}\n")
                fp.write("V:\n")
                if J == 0:
                    V = np.array([self.potential_local(S, S, S, J, Tz, ri) for ri in MeshPoints])
                    np.savetxt(fp, V, fmt="%.17f")
                else:
                    if not coup:
                        V = np.array([self.potential_local(J, J, S, J, Tz, ri) for ri in MeshPoints])
                        np.savetxt(fp, V, fmt="%.17f")
                    else:
                        Vpp = np.array([self.potential_local(J + 1, J + 1, S, J, Tz, ri) for ri in MeshPoints])
                        Vpm = np.array([self.potential_local(J + 1, J - 1, S, J, Tz, ri) for ri in MeshPoints])
                        Vmp = np.array([self.potential_local(J - 1, J + 1, S, J, Tz, ri) for ri in MeshPoints])
                        Vmm = np.array([self.potential_local(J - 1, J - 1, S, J, Tz, ri) for ri in MeshPoints])
                        V = np.block([[Vmm, Vmp], [Vpm, Vpp]])
                        np.savetxt(fp, V, fmt="%.17f")
        binary_file = "./" + dir + "/" + f"{self.chiral_type}_rmax{rmax}_N{N}_Jmax{Jmax}.bin"
        t2 = time.time()
        profiler.add_timing("Writing .dat", t2 - t1)
        # self.convert_dat_to_binary(file_name, binary_file)
        t3 = time.time()
        profiler.add_timing("Writing .bin", t3 - t2)
        print("files are generated to :\n", file_name, "\n", binary_file)

    @staticmethod
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
