## a realization of chiral NN potential under plane-wave basis,
## which is the input of nuclear matter calculation.


import numpy as np
import time


gA = 1.29
fpi = 92.4
Mp = 938.2720
Mn = 939.5654
Mnuc = 938.9183
Mpi = 138.0390
Mpi_charged = 139.5702
Mpi_neutral = 134.9766
hbarc = 197.32698

twopicubic = 248.0502134423986

lecs = {}
lecs["Lambda"] = 500
lecs["C_S_nn"] = -0.010011
lecs["C_T_nn"] = 0.000543


# automized spin projection.
def potential_auto(
    f1: float,
    f2: float,
    f3: float,
    f4: float,
    f5: float,
    f6: float,
    ppx: float,
    ppy: float,
    ppz: float,
    px: float,
    py: float,
    pz: float,
    s1_out: int,
    s2_out: int,
    s1_in: int,
    s2_in: int,
):
    imunit = 0
    spin_tuple = (s1_out, s2_out, s1_in, s2_in)
    if spin_tuple == (1, 1, 1, 1):
        return (
            f1
            + f2
            + f5 * np.power(ppz, 2)
            + f6 * np.power(ppz, 2)
            + f4 * np.power(ppy, 2) * np.power(px, 2)
            - 2 * f4 * ppx * ppy * px * py
            + f4 * np.power(ppx, 2) * np.power(py, 2)
            + imunit * (2 * f3 * ppy * px - 2 * f3 * ppx * py)
            + 2 * f5 * ppz * pz
            - 2 * f6 * ppz * pz
            + f5 * np.power(pz, 2)
            + f6 * np.power(pz, 2)
        )
    elif spin_tuple == (1, 1, 1, -1):
        return (
            f5 * ppx * ppz
            + f6 * ppx * ppz
            - f3 * ppz * px
            + f5 * ppz * px
            - f6 * ppz * px
            + f4 * ppy * ppz * px * py
            - f4 * ppx * ppz * np.power(py, 2)
            + f3 * ppx * pz
            + f5 * ppx * pz
            - f6 * ppx * pz
            + f5 * px * pz
            + f6 * px * pz
            - f4 * np.power(ppy, 2) * px * pz
            + f4 * ppx * ppy * py * pz
            + imunit
            * (
                -(f5 * ppy * ppz)
                - f6 * ppy * ppz
                + f4 * ppy * ppz * np.power(px, 2)
                + f3 * ppz * py
                - f5 * ppz * py
                + f6 * ppz * py
                - f4 * ppx * ppz * px * py
                - f3 * ppy * pz
                - f5 * ppy * pz
                + f6 * ppy * pz
                - f4 * ppx * ppy * px * pz
                - f5 * py * pz
                - f6 * py * pz
                + f4 * np.power(ppx, 2) * py * pz
            )
        )
    elif spin_tuple == (1, 1, -1, 1):
        return (
            f5 * ppx * ppz
            + f6 * ppx * ppz
            - f3 * ppz * px
            + f5 * ppz * px
            - f6 * ppz * px
            + f4 * ppy * ppz * px * py
            - f4 * ppx * ppz * np.power(py, 2)
            + f3 * ppx * pz
            + f5 * ppx * pz
            - f6 * ppx * pz
            + f5 * px * pz
            + f6 * px * pz
            - f4 * np.power(ppy, 2) * px * pz
            + f4 * ppx * ppy * py * pz
            + imunit
            * (
                -(f5 * ppy * ppz)
                - f6 * ppy * ppz
                + f4 * ppy * ppz * np.power(px, 2)
                + f3 * ppz * py
                - f5 * ppz * py
                + f6 * ppz * py
                - f4 * ppx * ppz * px * py
                - f3 * ppy * pz
                - f5 * ppy * pz
                + f6 * ppy * pz
                - f4 * ppx * ppy * px * pz
                - f5 * py * pz
                - f6 * py * pz
                + f4 * np.power(ppx, 2) * py * pz
            )
        )
    elif spin_tuple == (1, 1, -1, -1):
        return (
            f5 * np.power(ppx, 2)
            + f6 * np.power(ppx, 2)
            - f5 * np.power(ppy, 2)
            - f6 * np.power(ppy, 2)
            + 2 * f5 * ppx * px
            - 2 * f6 * ppx * px
            + f5 * np.power(px, 2)
            + f6 * np.power(px, 2)
            - f4 * np.power(ppz, 2) * np.power(px, 2)
            - 2 * f5 * ppy * py
            + 2 * f6 * ppy * py
            - f5 * np.power(py, 2)
            - f6 * np.power(py, 2)
            + f4 * np.power(ppz, 2) * np.power(py, 2)
            + 2 * f4 * ppx * ppz * px * pz
            - 2 * f4 * ppy * ppz * py * pz
            - f4 * np.power(ppx, 2) * np.power(pz, 2)
            + f4 * np.power(ppy, 2) * np.power(pz, 2)
            + imunit
            * (
                -2 * f5 * ppx * ppy
                - 2 * f6 * ppx * ppy
                - 2 * f5 * ppy * px
                + 2 * f6 * ppy * px
                - 2 * f5 * ppx * py
                + 2 * f6 * ppx * py
                - 2 * f5 * px * py
                - 2 * f6 * px * py
                + 2 * f4 * np.power(ppz, 2) * px * py
                - 2 * f4 * ppy * ppz * px * pz
                - 2 * f4 * ppx * ppz * py * pz
                + 2 * f4 * ppx * ppy * np.power(pz, 2)
            )
        )
    elif spin_tuple == (1, -1, 1, 1):
        return (
            f5 * ppx * ppz
            + f6 * ppx * ppz
            + f3 * ppz * px
            + f5 * ppz * px
            - f6 * ppz * px
            + f4 * ppy * ppz * px * py
            - f4 * ppx * ppz * np.power(py, 2)
            - f3 * ppx * pz
            + f5 * ppx * pz
            - f6 * ppx * pz
            + f5 * px * pz
            + f6 * px * pz
            - f4 * np.power(ppy, 2) * px * pz
            + f4 * ppx * ppy * py * pz
            + imunit
            * (
                f5 * ppy * ppz
                + f6 * ppy * ppz
                - f4 * ppy * ppz * np.power(px, 2)
                + f3 * ppz * py
                + f5 * ppz * py
                - f6 * ppz * py
                + f4 * ppx * ppz * px * py
                - f3 * ppy * pz
                + f5 * ppy * pz
                - f6 * ppy * pz
                + f4 * ppx * ppy * px * pz
                + f5 * py * pz
                + f6 * py * pz
                - f4 * np.power(ppx, 2) * py * pz
            )
        )
    elif spin_tuple == (1, -1, 1, -1):
        return (
            f1
            - f2
            - f5 * np.power(ppz, 2)
            - f6 * np.power(ppz, 2)
            - f4 * np.power(ppy, 2) * np.power(px, 2)
            + 2 * f4 * ppx * ppy * px * py
            - f4 * np.power(ppx, 2) * np.power(py, 2)
            - 2 * f5 * ppz * pz
            + 2 * f6 * ppz * pz
            - f5 * np.power(pz, 2)
            - f6 * np.power(pz, 2)
        )
    elif spin_tuple == (1, -1, -1, 1):
        return (
            2 * f2
            + f5 * np.power(ppx, 2)
            + f6 * np.power(ppx, 2)
            + f5 * np.power(ppy, 2)
            + f6 * np.power(ppy, 2)
            + 2 * f5 * ppx * px
            - 2 * f6 * ppx * px
            + f5 * np.power(px, 2)
            + f6 * np.power(px, 2)
            + f4 * np.power(ppz, 2) * np.power(px, 2)
            + 2 * f5 * ppy * py
            - 2 * f6 * ppy * py
            + f5 * np.power(py, 2)
            + f6 * np.power(py, 2)
            + f4 * np.power(ppz, 2) * np.power(py, 2)
            - 2 * f4 * ppx * ppz * px * pz
            - 2 * f4 * ppy * ppz * py * pz
            + f4 * np.power(ppx, 2) * np.power(pz, 2)
            + f4 * np.power(ppy, 2) * np.power(pz, 2)
        )
    elif spin_tuple == (1, -1, -1, -1):
        return (
            -(f5 * ppx * ppz)
            - f6 * ppx * ppz
            - f3 * ppz * px
            - f5 * ppz * px
            + f6 * ppz * px
            - f4 * ppy * ppz * px * py
            + f4 * ppx * ppz * np.power(py, 2)
            + f3 * ppx * pz
            - f5 * ppx * pz
            + f6 * ppx * pz
            - f5 * px * pz
            - f6 * px * pz
            + f4 * np.power(ppy, 2) * px * pz
            - f4 * ppx * ppy * py * pz
            + imunit
            * (
                f5 * ppy * ppz
                + f6 * ppy * ppz
                - f4 * ppy * ppz * np.power(px, 2)
                + f3 * ppz * py
                + f5 * ppz * py
                - f6 * ppz * py
                + f4 * ppx * ppz * px * py
                - f3 * ppy * pz
                + f5 * ppy * pz
                - f6 * ppy * pz
                + f4 * ppx * ppy * px * pz
                + f5 * py * pz
                + f6 * py * pz
                - f4 * np.power(ppx, 2) * py * pz
            )
        )
    elif spin_tuple == (-1, 1, 1, 1):
        return (
            f5 * ppx * ppz
            + f6 * ppx * ppz
            + f3 * ppz * px
            + f5 * ppz * px
            - f6 * ppz * px
            + f4 * ppy * ppz * px * py
            - f4 * ppx * ppz * np.power(py, 2)
            - f3 * ppx * pz
            + f5 * ppx * pz
            - f6 * ppx * pz
            + f5 * px * pz
            + f6 * px * pz
            - f4 * np.power(ppy, 2) * px * pz
            + f4 * ppx * ppy * py * pz
            + imunit
            * (
                f5 * ppy * ppz
                + f6 * ppy * ppz
                - f4 * ppy * ppz * np.power(px, 2)
                + f3 * ppz * py
                + f5 * ppz * py
                - f6 * ppz * py
                + f4 * ppx * ppz * px * py
                - f3 * ppy * pz
                + f5 * ppy * pz
                - f6 * ppy * pz
                + f4 * ppx * ppy * px * pz
                + f5 * py * pz
                + f6 * py * pz
                - f4 * np.power(ppx, 2) * py * pz
            )
        )
    elif spin_tuple == (-1, 1, 1, -1):
        return (
            2 * f2
            + f5 * np.power(ppx, 2)
            + f6 * np.power(ppx, 2)
            + f5 * np.power(ppy, 2)
            + f6 * np.power(ppy, 2)
            + 2 * f5 * ppx * px
            - 2 * f6 * ppx * px
            + f5 * np.power(px, 2)
            + f6 * np.power(px, 2)
            + f4 * np.power(ppz, 2) * np.power(px, 2)
            + 2 * f5 * ppy * py
            - 2 * f6 * ppy * py
            + f5 * np.power(py, 2)
            + f6 * np.power(py, 2)
            + f4 * np.power(ppz, 2) * np.power(py, 2)
            - 2 * f4 * ppx * ppz * px * pz
            - 2 * f4 * ppy * ppz * py * pz
            + f4 * np.power(ppx, 2) * np.power(pz, 2)
            + f4 * np.power(ppy, 2) * np.power(pz, 2)
        )
    elif spin_tuple == (-1, 1, -1, 1):
        return (
            f1
            - f2
            - f5 * np.power(ppz, 2)
            - f6 * np.power(ppz, 2)
            - f4 * np.power(ppy, 2) * np.power(px, 2)
            + 2 * f4 * ppx * ppy * px * py
            - f4 * np.power(ppx, 2) * np.power(py, 2)
            - 2 * f5 * ppz * pz
            + 2 * f6 * ppz * pz
            - f5 * np.power(pz, 2)
            - f6 * np.power(pz, 2)
        )
    elif spin_tuple == (-1, 1, -1, -1):
        return (
            -(f5 * ppx * ppz)
            - f6 * ppx * ppz
            - f3 * ppz * px
            - f5 * ppz * px
            + f6 * ppz * px
            - f4 * ppy * ppz * px * py
            + f4 * ppx * ppz * np.power(py, 2)
            + f3 * ppx * pz
            - f5 * ppx * pz
            + f6 * ppx * pz
            - f5 * px * pz
            - f6 * px * pz
            + f4 * np.power(ppy, 2) * px * pz
            - f4 * ppx * ppy * py * pz
            + imunit
            * (
                f5 * ppy * ppz
                + f6 * ppy * ppz
                - f4 * ppy * ppz * np.power(px, 2)
                + f3 * ppz * py
                + f5 * ppz * py
                - f6 * ppz * py
                + f4 * ppx * ppz * px * py
                - f3 * ppy * pz
                + f5 * ppy * pz
                - f6 * ppy * pz
                + f4 * ppx * ppy * px * pz
                + f5 * py * pz
                + f6 * py * pz
                - f4 * np.power(ppx, 2) * py * pz
            )
        )
    elif spin_tuple == (-1, -1, 1, 1):
        return (
            f5 * np.power(ppx, 2)
            + f6 * np.power(ppx, 2)
            - f5 * np.power(ppy, 2)
            - f6 * np.power(ppy, 2)
            + 2 * f5 * ppx * px
            - 2 * f6 * ppx * px
            + f5 * np.power(px, 2)
            + f6 * np.power(px, 2)
            - f4 * np.power(ppz, 2) * np.power(px, 2)
            - 2 * f5 * ppy * py
            + 2 * f6 * ppy * py
            - f5 * np.power(py, 2)
            - f6 * np.power(py, 2)
            + f4 * np.power(ppz, 2) * np.power(py, 2)
            + 2 * f4 * ppx * ppz * px * pz
            - 2 * f4 * ppy * ppz * py * pz
            - f4 * np.power(ppx, 2) * np.power(pz, 2)
            + f4 * np.power(ppy, 2) * np.power(pz, 2)
            + imunit
            * (
                2 * f5 * ppx * ppy
                + 2 * f6 * ppx * ppy
                + 2 * f5 * ppy * px
                - 2 * f6 * ppy * px
                + 2 * f5 * ppx * py
                - 2 * f6 * ppx * py
                + 2 * f5 * px * py
                + 2 * f6 * px * py
                - 2 * f4 * np.power(ppz, 2) * px * py
                + 2 * f4 * ppy * ppz * px * pz
                + 2 * f4 * ppx * ppz * py * pz
                - 2 * f4 * ppx * ppy * np.power(pz, 2)
            )
        )
    elif spin_tuple == (-1, -1, 1, -1):
        return (
            -(f5 * ppx * ppz)
            - f6 * ppx * ppz
            + f3 * ppz * px
            - f5 * ppz * px
            + f6 * ppz * px
            - f4 * ppy * ppz * px * py
            + f4 * ppx * ppz * np.power(py, 2)
            - f3 * ppx * pz
            - f5 * ppx * pz
            + f6 * ppx * pz
            - f5 * px * pz
            - f6 * px * pz
            + f4 * np.power(ppy, 2) * px * pz
            - f4 * ppx * ppy * py * pz
            + imunit
            * (
                -(f5 * ppy * ppz)
                - f6 * ppy * ppz
                + f4 * ppy * ppz * np.power(px, 2)
                + f3 * ppz * py
                - f5 * ppz * py
                + f6 * ppz * py
                - f4 * ppx * ppz * px * py
                - f3 * ppy * pz
                - f5 * ppy * pz
                + f6 * ppy * pz
                - f4 * ppx * ppy * px * pz
                - f5 * py * pz
                - f6 * py * pz
                + f4 * np.power(ppx, 2) * py * pz
            )
        )
    elif spin_tuple == (-1, -1, -1, 1):
        return (
            -(f5 * ppx * ppz)
            - f6 * ppx * ppz
            + f3 * ppz * px
            - f5 * ppz * px
            + f6 * ppz * px
            - f4 * ppy * ppz * px * py
            + f4 * ppx * ppz * np.power(py, 2)
            - f3 * ppx * pz
            - f5 * ppx * pz
            + f6 * ppx * pz
            - f5 * px * pz
            - f6 * px * pz
            + f4 * np.power(ppy, 2) * px * pz
            - f4 * ppx * ppy * py * pz
            + imunit
            * (
                -(f5 * ppy * ppz)
                - f6 * ppy * ppz
                + f4 * ppy * ppz * np.power(px, 2)
                + f3 * ppz * py
                - f5 * ppz * py
                + f6 * ppz * py
                - f4 * ppx * ppz * px * py
                - f3 * ppy * pz
                - f5 * ppy * pz
                + f6 * ppy * pz
                - f4 * ppx * ppy * px * pz
                - f5 * py * pz
                - f6 * py * pz
                + f4 * np.power(ppx, 2) * py * pz
            )
        )
    elif spin_tuple == (-1, -1, -1, -1):
        return (
            f1
            + f2
            + f5 * np.power(ppz, 2)
            + f6 * np.power(ppz, 2)
            + f4 * np.power(ppy, 2) * np.power(px, 2)
            - 2 * f4 * ppx * ppy * px * py
            + f4 * np.power(ppx, 2) * np.power(py, 2)
            + imunit * (-2 * f3 * ppy * px + 2 * f3 * ppx * py)
            + 2 * f5 * ppz * pz
            - 2 * f6 * ppz * pz
            + f5 * np.power(pz, 2)
            + f6 * np.power(pz, 2)
        )
    else:
        raise ValueError("Invalid spin configuration in function calculate_potential!")


def potential_lo_contact(
    ppx: float, ppy: float, ppz: float, px: float, py: float, pz: float, lecs: dict
):
    f1, f2, f3, f4, f5, f6 = [0, 0, 0, 0, 0, 0]
    pmag = np.sqrt(px**2 + py**2 + pz**2)
    ppmag = np.sqrt(ppx**2 + ppy**2 + ppz**2)

    # 10^4 GeV^-2 - > 1e-2 MeV^-2
    f1 = f1 + lecs["C_S_nn"] * 1e-2
    f2 = f2 + lecs["C_T_nn"] * 1e-2

    Lambda = lecs["Lambda"]
    n = 3
    regulator = np.exp(-(pmag ** (2 * n) + ppmag ** (2 * n)) / (Lambda ** (2 * n)))

    global Mnuc
    e = np.sqrt(Mnuc**2 + pmag**2)
    ep = np.sqrt(Mnuc**2 + ppmag**2)
    min_rel = np.sqrt(Mnuc / e) * np.sqrt(Mnuc / ep)

    return regulator * min_rel * np.array([f1, f2, f3, f4, f5, f6])


def potential_one_pion_exchange_neutral(
    ppx: float, ppy: float, ppz: float, px: float, py: float, pz: float, lecs: dict
):
    f1, f2, f3, f4, f5, f6 = [0, 0, 0, 0, 0, 0]
    pmag = np.sqrt(px**2 + py**2 + pz**2)
    ppmag = np.sqrt(ppx**2 + ppy**2 + ppz**2)
    qx = ppx - px
    qy = ppy - py
    qz = ppz - pz
    q2 = qx**2 + qy**2 + qz**2

    global gA, fpi, Mpi_neutral
    prefactor = -1.0 * gA**2 / (4 * fpi**2)
    f6 = f6 + prefactor / (q2 + Mpi_neutral**2)

    Lambda = lecs["Lambda"]
    n = 4
    regulator = np.exp(-(pmag ** (2 * n) + ppmag ** (2 * n)) / (Lambda ** (2 * n)))

    global Mnuc
    e = np.sqrt(Mnuc**2 + pmag**2)
    ep = np.sqrt(Mnuc**2 + ppmag**2)
    min_rel = np.sqrt(Mnuc / e) * np.sqrt(Mnuc / ep)

    return regulator * min_rel * np.array([f1, f2, f3, f4, f5, f6])


def potential_chiral_nn(
    p_out: np.array,
    s1_out: int,
    s2_out: int,
    p_in: np.array,
    s1_in: int,
    s2_in: int,
    lecs: dict,
):
    """
    The nn chiral potential under plane wave basis.
    param p_out: relative out momentum (length-3 array in MeV)
    param p_in : relative in momentum
    param s_i : +1 or -1, spin projection of particle i
    return: value of potential in MeV^(-2)
    note: to get discretized potential in MeV, multiply dk^3 or (2Pi/Lbox)^3
    """

    ppx, ppy, ppz = p_out[0], p_out[1], p_out[2]
    px, py, pz = p_in[0], p_in[1], p_in[2]
    f1, f2, f3, f4, f5, f6 = [0, 0, 0, 0, 0, 0]
    component_vec = np.array([f1, f2, f3, f4, f5, f6])

    component_vec = component_vec + potential_lo_contact(
        ppx, ppy, ppz, px, py, pz, lecs
    )
    component_vec = component_vec + potential_one_pion_exchange_neutral(
        ppx, ppy, ppz, px, py, pz, lecs
    )

    global twopicubic
    f1, f2, f3, f4, f5, f6 = component_vec
    mtx3d = np.real(
        potential_auto(
            f1,
            f2,
            f3,
            f4,
            f5,
            f6,
            ppx,
            ppy,
            ppz,
            px,
            py,
            pz,
            s1_out,
            s2_out,
            s1_in,
            s2_in,
        )
        / twopicubic
    )
    return mtx3d


def main():
    global lecs
    p_out = np.array([100, 200, 300])
    p_in = np.array([0, 200, 300])
    mtx = potential_chiral_nn(p_out, -1, -1, p_in, -1, -1, lecs)
    print(mtx)
    # for s1_out in [-1, 1]:
    #     for s2_out in [-1, 1]:
    #         for s1_in in [-1, 1]:
    #             for s2_in in [-1, 1]:
    #                 # print(s1_out, s2_out, s1_in, s2_in)
    #                 mtx = potential_chiral_nn(
    #                     p_out, s1_out, s2_out, p_in, s1_in, s2_in, lecs
    #                 )
    #                 print(mtx)


if __name__ == "__main__":
    main()
