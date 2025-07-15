import sys

sys.path.append("./deps")

import math
import random
import numpy as np
from typing import List, Tuple, Iterable
import chiral_potential
import constants as const
import basic_math


# stochastic SRG flow equation for the two-nucleon potential, test in progress.
class sSRG:

    def __init__(self, params: dict):
        self.potential: chiral_potential = chiral_potential.two_nucleon_potential(params["potential_type"])
        self.num_mesh: int
        self.mesh_q: np.ndarray
        self.weight_q: np.ndarray
        self.setup_mesh(params["q_min"], params["q_max"], params["q_number"])
        self.Tkin: np.ndarray
        self.setup_kinetic()
        self.V0: np.ndarray
        self.setup_initial_potential()
        self.H0: np.ndarray = self.Tkin + self.V0
        self.E0: float = np.linalg.eigvalsh(self.H0).min()
        self.num_H0: float = np.sum(np.abs(self.H0))
        self.eta: np.ndarray
        self.setup_generator()

        self.is_initiator: bool = False
        self.target_walker_number: float = params["target_walker_number"]
        self.random_sampling: bool = params["random_sampling"]
        self.loops: int = params["loops"]
        self.d_tau: float = params["d_tau"]
        self.A: int = params["A"]
        self.xi: float = params["xi"]
        self.zeta: float = params["zeta"]
        self.initiator_threshold: int = params["initiator_threshold"]
        self.seed: int = params["seed"]
        self.min_spawn_num: float = 0
        self.min_walker_num: float = 0
        self.S: float = 0.0
        self.walkers: dict = {}  # key = (i, j); value = fij
        self.new_walkers: List[Tuple[int, int, bool, float]] = []
        self.tau_loops_trace: List[List[float]] = []
        self.number_loops_trace: List[List[float]] = []
        self.shift_loops_trace: List[List[float]] = []

    # commutator, arguments A and B matrices
    @staticmethod
    def commutator(A, B):
        return A @ B - B @ A

    # make the matrix symmetric
    @staticmethod
    def make_symmetric(A):
        return (A + A.T) / 2

    # q_min, q_max in fm^(-1)
    def setup_mesh(self, q_min: float, q_max: float, q_number: int):
        self.mesh_q = const.hbarc * np.linspace(q_min, q_max, q_number)
        self.weight_q = const.hbarc * np.full(q_number, (q_max - q_min) / q_number)
        self.num_mesh = q_number

    def setup_kinetic(self):
        Nq = self.num_mesh
        self.Tkin = np.zeros((Nq, Nq))
        for i in range(Nq):
            pi = self.mesh_q[i]
            self.Tkin[i, i] = pi**2 / const.MN
        print("finish setting up kinetic energy matrix.")

    def setup_initial_potential(self):
        # consider np 1S0 channel for now.
        ll, l, j, s, tz = 0, 0, 0, 0, 0
        Nq = self.num_mesh
        self.V0 = np.zeros((Nq, Nq))
        for idx_pp in range(Nq):
            for idx_p in range(Nq):
                pp = self.mesh_q[idx_pp]
                p = self.mesh_q[idx_p]
                wpp = self.weight_q[idx_pp]
                wp = self.weight_q[idx_p]
                self.V0[idx_pp, idx_p] = np.sqrt(wpp * wp) * pp * p * self.potential.potential(ll, l, pp, p, j, s, tz)
        print("finish setting up initial potential matrix.")

    def setup_generator(self):
        self.eta = self.commutator(self.Tkin, self.H0)
        print("finish setting up initial generator matrix.")

    # Dress the weights for mtx, MeV^(-2) to MeV
    def dress_weights(self, mtx: np.ndarray):
        Nq = self.num_mesh
        for idx_pp in range(Nq):
            for idx_p in range(Nq):
                pp = self.mesh_q[idx_pp]
                p = self.mesh_q[idx_p]
                wpp = self.weight_q[idx_pp]
                wp = self.weight_q[idx_p]
                mtx[idx_pp, idx_p] *= np.sqrt(wpp * wp) * pp * p
        return mtx

    # Undress the weights for mtx, MeV to MeV^(-2)
    def undress_weights(self, mtx: np.ndarray):
        Nq = self.num_mesh
        for idx_pp in range(Nq):
            for idx_p in range(Nq):
                pp = self.mesh_q[idx_pp]
                p = self.mesh_q[idx_p]
                wpp = self.weight_q[idx_pp]
                wp = self.weight_q[idx_p]
                mtx[idx_pp, idx_p] /= np.sqrt(wpp * wp) * pp * p
        return mtx

    def get_number(self) -> float:
        sum_walkers = sum(abs(value) for value in self.walkers.values())
        return sum_walkers

    # update the generator matrix eta: eta = [T, H]
    def update_eta(self):
        H_now = self.get_H()
        self.eta = self.commutator(self.Tkin, H_now)

    def get_H(self) -> np.ndarray:
        Nq = self.num_mesh
        H = np.zeros((Nq, Nq))
        for (i, j), fij in self.walkers.items():
            H[i, j] += fij
        return self.num_H0 / self.target_walker_number * H

    def get_V(self) -> np.ndarray:
        V_now = self.get_H() - self.Tkin
        return self.undress_weights(V_now)

    def get_E(self):
        H_now = self.get_H()
        E_now = np.linalg.eigvalsh(H_now).min()
        return E_now

    def abs_cut_to(self, num: float, target: float) -> float:
        abs_num = abs(num)
        if abs_num >= target:
            return num
        else:
            sign = 1.0 if num >= 0 else -1.0
            return sign * target * float(random.uniform(0, target) < abs_num)

    def walker_num_cut(self, num: float) -> float:
        abs_num = abs(num)
        if abs_num >= self.min_walker_num:
            return num
        else:
            sign = 1.0 if num >= 0 else -1.0
            return sign * self.min_walker_num * float(random.uniform(0, self.min_walker_num) < abs_num)

    def normalize_walker_number(self):
        total_number = self.get_number()
        fac = self.target_walker_number / total_number
        for (i, j), fij in self.walkers.items():
            self.walkers[(i, j)] = fij * fac

    def initialize_walkers(self):
        if self.random_sampling:
            self.initialize_random_walkers()
        else:
            self.initialize_exact_walkers()

    def initialize_exact_walkers(self):
        Nq = self.num_mesh
        for i in range(Nq):
            for j in range(Nq):
                fij = self.H0[i, j]
                self.walkers[(i, j)] = fij
        self.normalize_walker_number()
        # print("finish initializing exact walkers.")

    # using importance sampling to initialize walkers
    def initialize_random_walkers(self):
        Nq = self.num_mesh
        candidate_i = np.random.randint(0, Nq, self.target_walker_number)
        candidate_j = np.random.randint(0, Nq, self.target_walker_number)
        weights = self.H0[candidate_i, candidate_j]
        normalized_weights = weights / np.sum(weights)
        chosen_indices = np.random.choice(np.arange(self.target_walker_number), size=self.target_walker_number, p=normalized_weights, replace=True)
        walkers_sir_i = candidate_i[chosen_indices]
        walkers_sir_j = candidate_j[chosen_indices]
        for idx in range(self.target_walker_number):
            i = walkers_sir_i[idx]
            j = walkers_sir_j[idx]
            Hij = self.H0[i, j]
            if (i, j) not in self.walkers:
                self.walkers[(i, j)] = math.copysign(1.0, Hij)
            else:
                self.walkers[(i, j)] += math.copysign(1.0, Hij)
        # print("finish initializing random walkers.")

    def make_walkers_symmetric(self):
        for (i, j), cij in list(self.walkers.items()):
            if i == j:
                continue
            if (j, i) not in self.walkers:
                self.walkers[(i, j)] = cij / 2
                self.walkers[(j, i)] = cij / 2
            else:
                num_average = (self.walkers[(i, j)] + self.walkers[(j, i)]) / 2
                self.walkers[(i, j)] = num_average
                self.walkers[(j, i)] = num_average

    def print_walkers(self):
        print("walkers:")
        for (i, j), cij in self.walkers.items():
            print(f"({i}, {j}): {cij:.3e}")
        print("total walkers:", len(self.walkers))
        print("total walker number:", self.get_number())

    def step(self):
        self.seed += 1
        random.seed(self.seed)
        for (i, j), cij in list(self.walkers.items()):
            Cij = self.walker_num_cut(cij)
            if Cij == 0.0:
                del self.walkers[(i, j)]
                continue
            Nij = math.floor(Cij + random.uniform(0, 1))
            eta_ii = self.eta[i, i]
            eta_jj = self.eta[j, j]
            # pd = self.d_tau * (-eta_ii + eta_jj - self.S)
            pd = self.d_tau * (-eta_ii + eta_jj)
            self.walkers[(i, j)] = cij - pd * Nij  # diagonal step
            # is_initiator = self.is_initiator and (abs(Cij) > self.initiator_threshold)
            is_initiator = True
            for _ in range(abs(Nij)):
                k, l = 0, 0
                while True:
                    k = random.randint(0, self.num_mesh - 1)
                    l = random.randint(0, self.num_mesh - 1)
                    if (k, l) != (i, j):
                        break
                eta_ki = self.eta[k, i]
                eta_lj = self.eta[l, j]
                invp_ik = self.num_mesh - 1
                invp_jl = self.num_mesh - 1
                spawn_num_ik = basic_math.sign(Nij) * self.d_tau * eta_ki * invp_ik
                spawn_num_jl = -basic_math.sign(Nij) * self.d_tau * eta_lj * invp_jl
                # spawn_num_ik = self.abs_cut_to(spawn_num_ik, self.min_spawn_num)
                # spawn_num_jl = self.abs_cut_to(spawn_num_jl, self.min_spawn_num)
                if spawn_num_ik == 0.0 and spawn_num_jl == 0.0:
                    continue
                # spawning step
                self.new_walkers.append((k, j, is_initiator, spawn_num_ik))
                self.new_walkers.append((i, l, is_initiator, spawn_num_jl))

    def annihilation(self):
        self.seed += 1
        random.seed(self.seed)
        for i, j, is_initiator, spawn_num in self.new_walkers:
            if ((i, j) not in self.walkers) and (not is_initiator):
                continue
            else:
                if (i, j) in self.walkers:
                    current_num = self.walkers[(i, j)]
                    new_num = current_num + spawn_num
                    if new_num == 0.0:
                        del self.walkers[(i, j)]
                        continue
                    else:
                        self.walkers[(i, j)] = new_num
                else:
                    new_num = spawn_num
                    self.walkers[(i, j)] = new_num
        self.new_walkers.clear()

    def start(self, tau_target: float):
        steps = int(tau_target / self.d_tau)
        print(f"total steps per loop = {steps}")
        self.d_tau = tau_target / float(steps)
        if steps < 1:
            raise ValueError("error: steps < 1 in sSRG!")
        print("! evolution begins")
        print(f"!{'loop':>5}{'step':>12}{'S':>16}{'E':>16}{'Nw':>16}")
        for l in range(self.loops):
            self.walkers.clear()
            self.initialize_walkers()
            new_num = self.get_number()
            old_num = new_num
            tau_trace_this_loop = []
            number_trace_this_loop = []
            shift_trace_this_loop = []
            for i in range(1, steps + 1):
                tau_now = self.d_tau * i
                self.step()
                self.annihilation()
                self.make_walkers_symmetric()
                self.update_eta()
                if i % self.A == 0:
                    old_num = new_num
                    new_num = self.get_number()
                    energy_now = self.get_E()
                    print(f"{l:>6}{i:>12}{self.S:>16.3f}{energy_now:>16.3e}{new_num:>16.3e}")
                    # self.S = self.S - self.xi / (self.A * self.d_tau) * math.log(new_num / old_num) - self.zeta / (self.A * self.d_tau) * math.log(new_num / self.target_walker_number)
                    tau_trace_this_loop.append(tau_now)
                    number_trace_this_loop.append(new_num)
                    shift_trace_this_loop.append(self.S)
            self.tau_loops_trace.append(tau_trace_this_loop)
            self.number_loops_trace.append(number_trace_this_loop)
            self.shift_loops_trace.append(shift_trace_this_loop)
            print("! evolution ends")
