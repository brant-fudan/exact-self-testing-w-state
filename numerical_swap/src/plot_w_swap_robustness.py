"""Publication figures from audited W-state SWAP results; no fitted curves."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


STYLE = {
    "marginal": {"color": "#0072B2", "marker": "o", "label": r"Marginal $Z,X,D$ (216 words)"},
    "joint": {"color": "#D55E00", "marker": "s", "label": r"Joint $[2,2,2]$ (125 words)"},
}


def load_rows(root):
    rows, inputs = [], {}
    for protocol in ("marginal", "joint"):
        for model in ("box", "white"):
            path = root / f"{protocol}_{model}" / "curve.csv"
            inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            with path.open() as handle:
                part = list(csv.DictReader(handle))
            if not part or any(r["diagnostics_pass"] != "True" for r in part):
                raise ValueError(f"Missing points or failed diagnostics: {path}")
            for row in part:
                for key in ("epsilon", "numerical_lower_bound"):
                    row[key] = float(row[key])
            rows.extend(part)
    return rows, inputs


def panel(ax, rows, model, full=False):
    for protocol, style in STYLE.items():
        data = sorted((r for r in rows if r["protocol"] == protocol and r["noise_model"] == model
                       and (full or r["epsilon"] <= .01)), key=lambda r: r["epsilon"])
        ax.plot([r["epsilon"] for r in data], [r["numerical_lower_bound"] for r in data],
                linewidth=1.35, markersize=3.8, markeredgewidth=.65,
                markerfacecolor="white", **style)
    ax.axhline(2 / 3, color="#737373", ls=(0, (4, 3)), lw=.9, zorder=0)
    if full:
        ax.set_xscale("log")
        ax.set_xlim(8e-6, .065)
    else:
        ax.scatter([0], [1], marker="*", s=64, facecolors="white", edgecolors="black",
                   linewidths=.8, zorder=5, label="Exact ideal theorem")
        ax.set_xlim(-.00025, .01025)
        ax.set_xticks(np.linspace(0, .01, 6))
        ax.set_xticklabels(["0", "0.002", "0.004", "0.006", "0.008", "0.010"])
    ax.set_ylim(-.025, 1.025)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_xlabel(r"Absolute correlator deviation $\varepsilon$" if model == "box"
                  else r"White-noise parameter $p$")
    ax.grid(axis="y", alpha=.17, linewidth=.55)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", width=.7, length=3)


def save(fig, base):
    for extension in ("pdf", "svg", "png"):
        fig.savefig(base.with_suffix("." + extension), dpi=450, facecolor="white",
                    bbox_inches="tight", pad_inches=.04)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    rows, hashes = load_rows(args.results)
    output = args.results / "figures"
    output.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.labelsize": 8.5, "axes.titlesize": 9,
                         "axes.linewidth": .75, "legend.fontsize": 7.3,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none", "mathtext.fontset": "dejavusans"})
    for full, name in ((False, "w_swap_robustness"), (True, "w_swap_robustness_full_range")):
        fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.05), sharey=True)
        for ax, model, title in zip(axes, ("box", "white"),
                                     ("(a) Absolute-error box", "(b) White-noise path")):
            panel(ax, rows, model, full)
            ax.set_title(title, loc="left", pad=9)
        axes[0].set_ylabel("Numerical SWAP fidelity lower bound")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.52, .99),
                   frameon=False, ncol=3 if not full else 2)
        fig.text(.52, .018, r"Dashed line: biseparable threshold $F=2/3$. Lines only guide the eye.",
                 ha="center", fontsize=7.2, color="#444444")
        fig.subplots_adjust(left=.085, right=.99, bottom=.19, top=.79, wspace=.16)
        save(fig, output / name)
    for model, name in (("box", "w_swap_absolute_deviation"), ("white", "w_swap_white_noise")):
        fig, ax = plt.subplots(figsize=(3.6, 2.8))
        panel(ax, rows, model)
        ax.set_ylabel("Numerical SWAP fidelity lower bound")
        ax.legend(loc="lower left", frameon=False, fontsize=7)
        ax.text(.0099, 2 / 3 + .018, r"$F=2/3$", ha="right", fontsize=7, color="#666666")
        fig.tight_layout(pad=.6)
        save(fig, output / name)
    with (args.results / "all_curves.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "provenance.json").write_text(json.dumps({"input_sha256": hashes,
        "plot_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "value_column": "numerical_lower_bound", "interpolation": "none; straight guides only",
        "ideal_endpoint": "exact theorem, separate star; not overwritten numerical data",
        "diagnostics_required": True}, indent=2) + "\n")
    print(f"Saved 4 figure sets (PDF/SVG/450-dpi PNG) and {len(rows)} audited points in {args.results}")


if __name__ == "__main__":
    main()
