#!/usr/bin/env python3
"""
Run rust_ising_example for each grid size and plot CDFs.

Usage:
    python3 run_ising_bench.py --only 2x2,3x3        # test small sizes
    python3 run_ising_bench.py --skip 2x2,3x3,4x4    # skip already-done
    python3 run_ising_bench.py                        # run all
"""
import argparse
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/home/seb/rust_ising_and_graphical_tools/graphical_tools")
BINARY = ROOT / "testing/rust_ising_example"
INSTANCES_DIR = ROOT / "testing/ising_instances"
OUTPUT_DIR = ROOT / "testing/results"
PLOTTER = ROOT / "plot_cdf.py"

B_CAP = 10_000_000_000

JOBS = [
    # (grid_label, sizes_tuple, runs, target_fraction)
    ("2x2", (2, 2), 1000, 0.75),
    ("3x3", (3, 3), 1000, 0.75),
    ("4x4", (4, 4), 1000, 0.75),
    ("5x5", (5, 5), 1000, 0.75),
    ("6x5", (6, 5), 100, 0.10),
]


def compute_max_steps(N: int, target: float) -> int:
    p = 2.0 ** (-N)
    if p < 1e-16:
        t = -math.log(1 - target) * (2**N)
    else:
        t = math.log(1 - target) / math.log(1 - p)
    rounded = max(10, round(t / 10) * 10)
    return min(rounded, B_CAP)


def run_one(label: str, sizes: tuple[int, ...], runs: int, target: float, dry_run: bool = False) -> bool:
    N = sizes[0] * sizes[1]
    max_steps = compute_max_steps(N, target)
    input_file = INSTANCES_DIR / f"ising_{label}.txt"
    hits_file = OUTPUT_DIR / f"hits_{label}.txt"
    pdf_file = OUTPUT_DIR / f"cdf_{label}.pdf"
    csv_file = OUTPUT_DIR / f"stats_{label}.csv"

    if not input_file.exists():
        print(f"  SKIP: missing input {input_file}", file=sys.stderr)
        return False

    print(f"\n{'=' * 60}")
    print(f"  {label}: N={N}, max_steps={max_steps:,}, runs={runs}")
    print(f"{'=' * 60}")

    if dry_run:
        print(f"  [DRY RUN] would run rust_ising then plot")
        return True

    # Phase 1: run rust_ising
    t0 = time.time()
    result = subprocess.run(
        [str(BINARY),
         "-n", str(max_steps),
         "-r", str(runs),
         "-s", str(hits_file),
         str(input_file)],
        capture_output=True, text=True,
    )
    elapsed = time.time() - t0

    if result.returncode != 0:
        print(f"  FAILED ({elapsed:.1f}s): {result.stderr}", file=sys.stderr)
        return False

    print(f"  rust_ising done in {elapsed:.1f}s")

    # Phase 2: plot CDF
    result2 = subprocess.run(
        ["python3", str(PLOTTER),
         str(hits_file),
         "-N", str(N),
         "-o", str(pdf_file),
         "--csv", str(csv_file),
         "--style", "dark",
         "--max-steps", str(max_steps),
         "--title", f"Ising {label} ({N} spins) — Hitting Time CDF"],
        capture_output=True, text=True,
    )
    print(f"  plot: {result2.stdout.strip()}")
    if result2.returncode != 0:
        print(f"  plot FAILED: {result2.stderr}", file=sys.stderr)
        return False

    return True


def main():
    parser = argparse.ArgumentParser(description="Run Ising benchmarks")
    parser.add_argument("--only", type=str, default=None,
                        help="Comma-separated labels to run (e.g. 2x2,3x3)")
    parser.add_argument("--skip", type=str, default=None,
                        help="Comma-separated labels to skip")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would run without executing")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(exist_ok=True)

    only = set(x.strip() for x in (args.only or "").split(",") if x.strip()) if args.only else None
    skip = set(x.strip() for x in (args.skip or "").split(",") if x.strip()) if args.skip else set()

    jobs = JOBS
    if only:
        jobs = [j for j in jobs if j[0] in only]

    print(f"Running {len(jobs)} job(s), saving to {OUTPUT_DIR}/")
    total_ok = 0
    total_fail = 0

    for label, sizes, runs, target in jobs:
        if label in skip:
            print(f"\n  {label}: SKIPPED")
            continue
        ok = run_one(label, sizes, runs, target, dry_run=args.dry_run)
        if ok:
            total_ok += 1
        else:
            total_fail += 1

    print(f"\nDone: {total_ok} ok, {total_fail} failed")


if __name__ == "__main__":
    main()