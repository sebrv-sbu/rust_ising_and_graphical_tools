#!/usr/bin/env python3
"""
Plot empirical CDF of hitting times overlaid with theoretical curve.

Reads a stationary-hits file and plots the empirical CDF alongside
the theoretical curve:  1 - (1 - 2^{-N})^t

Usage:
    python3 plot_cdf.py stationary_hits.txt -N 16
    python3 plot_cdf.py stationary_hits.txt -N 16 -o cdf.png --csv stats.csv
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def load_hitting_times(path: Path, max_steps: int | None = None) -> tuple[np.ndarray, int, int]:
    """Parse the stationary-hits file."""
    times = []
    na_count = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.upper() == "NA":
                na_count += 1
            else:
                try:
                    times.append(int(line))
                except ValueError:
                    print(f"Warning: skipping unparseable line: {line!r}", file=sys.stderr)

    times = np.array(sorted(times), dtype=int)
    n_runs = len(times) + na_count

    if max_steps is None:
        max_val = times[-1] if len(times) > 0 else 0
        max_steps = max_val

    return times, n_runs, max_steps


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot empirical and theoretical CDF of hitting times."
    )
    parser.add_argument(
        "input_file", type=Path,
        help="Stationary hits file (output of rust_ising)"
    )
    parser.add_argument(
        "-N", "--n-vertices", type=int, required=True,
        help="Number of vertices N for theoretical curve 1−(1−2^{−N})^t"
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Output path (default: show interactively); saved as PDF unless name ends in .png"
    )
    parser.add_argument(
        "--csv", type=Path, default=None,
        help="Output CSV path for stats (default: stdout)"
    )
    parser.add_argument(
        "--title", type=str, default="Ising Model — Hitting Time CDF",
        help="Plot title"
    )
    parser.add_argument(
        "--style", type=str, default="default",
        choices=["default", "dark"],
        help="Visual style (default or dark)"
    )
    parser.add_argument(
        "--dpi", type=int, default=150,
        help="Output image DPI"
    )
    parser.add_argument(
        "--max-steps", type=int, default=None,
        help="Simulation max_steps; x-axis cut at next multiple of 10 (default: power of 10 from data)"
    )
    args = parser.parse_args()

    hits, n_runs, max_observed = load_hitting_times(args.input_file)

    if n_runs == 0:
        print("Error: no runs found in input file.", file=sys.stderr)
        sys.exit(1)

    n_hits = len(hits)
    n_na = n_runs - n_hits

    # --- Empirical CDF ---
    unique_times, counts = np.unique(hits, return_counts=True)
    cum_counts = np.cumsum(counts)
    cdf_values = cum_counts / n_runs

    x_vals = np.concatenate([[0], unique_times])
    y_vals = np.concatenate([[0.0], cdf_values])

    # --- Theoretical curve ---
    N = args.n_vertices
    p_hit = 2.0 ** (-N)
    mean_theory = 1.0 / p_hit
    median_theory = np.log(0.5) / np.log(1.0 - p_hit)

    end_empirical = x_vals[-1] if len(x_vals) > 0 else 1
    if args.max_steps is not None:
        xlim_right = ((args.max_steps + 9) // 10) * 10  # next multiple of 10
    else:
        xlim_right = 10 ** np.ceil(np.log10(max(end_empirical, 1)))

    t_theory = np.logspace(0, np.log10(xlim_right), 500)
    cdf_theory = 1.0 - (1.0 - p_hit) ** t_theory

    # --- Style ---
    if args.style == "dark":
        plt.style.use("dark_background")
    else:
        plt.style.use("seaborn-v0_8-whitegrid")

    fig, ax = plt.subplots(figsize=(10, 5))

    # Empirical CDF
    ax.step(x_vals, y_vals, where="post", color="#2a7ae2", linewidth=2.2,
            label=f"Empirical CDF (n={n_runs})")
    ax.fill_between(x_vals, y_vals, step="post", alpha=0.15, color="#2a7ae2")

    # Theoretical curve
    ax.plot(t_theory, cdf_theory, color="#ff7f0e", linewidth=2.0, linestyle="--",
            label=f"$1-(1-2^{{-{N}}})^t$", zorder=10)

    # --- Axes ---
    ax.set_xlim(left=0, right=xlim_right * 1.01)
    ax.set_ylim(bottom=0, top=1.02)

    ax.set_xlabel("Step (t)", fontsize=12, fontweight="medium")
    ax.set_ylabel("P(hit ≤ step)", fontsize=12, fontweight="medium")
    ax.set_title(args.title, fontsize=14, fontweight="bold", pad=12)

    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5))

    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)

    plt.tight_layout()

    if args.output:
        if args.output.suffix == '.png':
            out = args.output
            fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
        else:
            out = args.output if args.output.suffix == '.pdf' else args.output.with_suffix('.pdf')
            fig.savefig(out, bbox_inches="tight")
        print(f"Saved to {out}")
    else:
        plt.show()

    # --- CSV output ---
    rows = [
        {"stat": "N", "value": N},
        {"stat": "p_hit", "value": p_hit},
        {"stat": "median_theoretical", "value": f"{median_theory:.2f}"},
        {"stat": "mean_theoretical", "value": f"{mean_theory:.2f}"},
        {"stat": "runs", "value": n_runs},
        {"stat": "hits", "value": n_hits},
        {"stat": "censored_na", "value": n_na},
    ]
    if n_hits > 0:
        rows += [
            {"stat": "median_empirical", "value": f"{np.median(hits):.2f}"},
            {"stat": "mean_empirical", "value": f"{hits.mean():.2f}"},
            {"stat": "min_empirical", "value": hits.min()},
            {"stat": "max_empirical", "value": hits.max()},
        ]

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["stat", "value"])
            writer.writeheader()
            writer.writerows(rows)
        print(f"Stats written to {args.csv}")
    else:
        writer = csv.DictWriter(sys.stdout, fieldnames=["stat", "value"])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()