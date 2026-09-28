import argparse
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# DEFAULT_INPUT = "result/srg-stoch-phase-1s0-n3loemn500-Lambda2.0-loop10-step1000-Nw10000.txt"
# DEFAULT_INPUT = "result/srg-stoch-phase-3p0-n3loemn500-Lambda2.0-loop10-step1000-Nw10000.txt"
DEFAULT_INPUT = "result/srg-stoch-phase-3sd1-n3loemn500-Lambda2.0-loop10-step10000-Nw100.txt"

def read_header(file_path):
    with open(file_path, "r") as f:
        return f.readline().lstrip("#").strip().split()


def default_output_path(input_path):
    path = Path(input_path)
    return path.with_name(path.stem + ".png")


def infer_exact_path(input_path):
    exact_path = Path(str(input_path).replace("srg-stoch-phase", "srg-exact-phase", 1))
    return exact_path if exact_path.exists() else None


def load_data(file_path):
    data = np.loadtxt(file_path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def plot_single_channel(data, label, exact_data=None):
    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    Tlab = data[:, 0]
    phase_mean = data[:, 1]
    phase_std = data[:, 2]

    ax.errorbar(
        Tlab,
        phase_mean,
        yerr=phase_std,
        fmt="o-",
        color="C0",
        ecolor="C0",
        capsize=3,
        markersize=4,
        linewidth=1.2,
        elinewidth=0.9,
        label="stochastic SRG",
    )
    if exact_data is not None:
        ax.plot(
            exact_data[:, 0],
            exact_data[:, 1],
            "s--",
            color="black",
            markersize=3.5,
            linewidth=1.2,
            label="exact SRG",
        )
    ax.set_xlabel(r"$T_\mathrm{lab}$ [MeV]")
    ax.set_ylabel(r"$\delta$ [deg]")
    ax.set_title(f"SRG phase shift comparison: {label}")
    ax.minorticks_on()
    ax.tick_params(direction="in", top=True, right=True)
    ax.tick_params(which="minor", direction="in", top=True, right=True)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def plot_coupled_channel(data, label, exact_data=None):
    names = [r"$\bar{\delta}_{J-1}$", r"$\bar{\delta}_{J+1}$", r"$\bar{\epsilon}_J$"]
    fig, axs = plt.subplots(1, 3, figsize=(9.0, 3.0), sharex=True)
    Tlab = data[:, 0]
    means = data[:, 1:4]
    stds = data[:, 4:7]

    for i, ax in enumerate(axs):
        ax.errorbar(
            Tlab,
            means[:, i],
            yerr=stds[:, i],
            fmt="o-",
            color=f"C{i}",
            ecolor=f"C{i}",
            capsize=3,
            markersize=4,
            linewidth=1.2,
            elinewidth=0.9,
            label="stochastic SRG",
        )
        if exact_data is not None:
            ax.plot(
                exact_data[:, 0],
                exact_data[:, i + 1],
                "s--",
                color="black",
                markersize=3.5,
                linewidth=1.2,
                label="exact SRG",
            )
        ax.set_title(names[i])
        ax.set_xlabel(r"$T_\mathrm{lab}$ [MeV]")
        ax.minorticks_on()
        ax.tick_params(direction="in", top=True, right=True)
        ax.tick_params(which="minor", direction="in", top=True, right=True)
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    axs[0].set_ylabel("phase [deg]")
    fig.suptitle(f"SRG phase shift comparison: {label}")
    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description="Compare stochastic SRG phase shifts with exact SRG results.")
    parser.add_argument("input", nargs="?", default=DEFAULT_INPUT, help="phase txt file saved by renorm_stochastic_*.py")
    parser.add_argument("--exact", default=None, help="exact-SRG phase txt file; default: infer from stochastic file name")
    parser.add_argument("--no-exact", action="store_true", help="plot stochastic result only")
    parser.add_argument("-o", "--output", default=None, help="output figure path")
    parser.add_argument("--show", action="store_true", help="show the figure after saving")
    parser.add_argument("--label", default=None, help="legend/title label")
    args = parser.parse_args()

    input_path = Path(args.input)
    data = load_data(input_path)
    header = read_header(input_path)
    label = args.label or input_path.stem
    output_path = Path(args.output) if args.output else default_output_path(input_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    exact_data = None
    if not args.no_exact:
        exact_path = Path(args.exact) if args.exact else infer_exact_path(input_path)
        if exact_path is not None and exact_path.exists():
            exact_data = load_data(exact_path)
            print(f"using exact SRG phase file: {exact_path}")
        else:
            print("exact SRG phase file not found; plotting stochastic result only")

    if data.shape[1] == 3 and header[:3] == ["Tlab", "phase_mean", "phase_std"]:
        fig = plot_single_channel(data, label, exact_data)
    elif data.shape[1] == 7:
        fig = plot_coupled_channel(data, label, exact_data)
    else:
        raise ValueError(f"Unsupported phase file format: columns={data.shape[1]}, header={header}")

    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    print(f"saved: {output_path}")
    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
