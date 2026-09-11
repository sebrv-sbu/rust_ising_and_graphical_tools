#!/usr/bin/env python3
"""
Generate an Ising model input file with a ground state of all zeros (all spins down).

The Rust Ising code represents spins as bool: false = spin-down (σ = −1), true = spin-up (σ = +1).
The Hamiltonian is H = −Σ J_ij σ_i σ_j − Σ h_i σ_i.

For all-false (all −1) to be the ground state, we need:
  - J_ij > 0  (ferromagnetic: aligned spins of either sign minimize energy)
  - h_i < 0   (external field points in the −1 direction, favouring false)

Usage:
    python3 generate_ising.py                    # defaults: 2D, 4x4, temp=1.0
    python3 generate_ising.py --dim 3 --sizes 4 6 8 --temp 0.5 --J 1.5 --h -0.3
    python3 generate_ising.py -o my_config.txt --seed 42
"""

import argparse
import random
import sys
from pathlib import Path


def generate_ising_file(
    dim: int,
    sizes: list[int],
    temp: float,
    J: float | list[float] | None,
    h: float | list[float] | None,
    mu: float,
    output_path: Path,
    seed: int | None,
) -> None:
    if seed is not None:
        random.seed(seed)

    n_points = 1
    for s in sizes:
        n_points *= s

    deg = 2 * dim  # each dimension contributes 2 neighbours (± direction)

    # --- Compute coordinates and neighbours (same as to_coord / neighbours in ising.rs) ---
    def to_coord(node: int) -> list[int]:
        coord = []
        n = node
        for s in sizes:
            coord.append(n % s)
            n //= s
        return coord

    def from_coord(coord: list[int]) -> int:
        base = 1
        node = 0
        for pos, s in zip(coord, sizes):
            node += pos * base
            base *= s
        return node

    def neighbours(node: int) -> list[int]:
        nbrs = []
        coord = to_coord(node)
        for i, s in enumerate(sizes):
            if s > 1:
                orig = coord[i]
                coord[i] = (orig + 1) % s
                nbrs.append(from_coord(coord))
                coord[i] = (orig + s - 1) % s
                nbrs.append(from_coord(coord))
                coord[i] = orig
            else:
                nbrs.append(node)
                nbrs.append(node)
        return nbrs

    # --- Resolve J and h to per-site or per-edge values ---
    if J is None:
        # Default: random positive couplings (ferromagnetic)
        J_weights = [random.uniform(0.5, 2.0) for _ in range(n_points)]
    elif isinstance(J, (int, float)):
        J_weights = [float(J)] * n_points
    else:
        J_weights = list(J)
        if len(J_weights) != n_points:
            print(
                f"Warning: J has {len(J_weights)} entries, expected {n_points}. "
                f"{'Truncating.' if len(J_weights) > n_points else 'Padding with last value.'}"
            )
            if len(J_weights) > n_points:
                J_weights = J_weights[:n_points]
            else:
                pad = J_weights[-1] if J_weights else 1.0
                J_weights.extend([pad] * (n_points - len(J_weights)))

    if h is None:
        # Default: negative field (favours spin-down / false)
        h_values = [random.uniform(-2.0, -0.1) for _ in range(n_points)]
    elif isinstance(h, (int, float)):
        h_values = [float(h)] * n_points
    else:
        h_values = list(h)
        if len(h_values) != n_points:
            print(
                f"Warning: h has {len(h_values)} entries, expected {n_points}. "
                f"{'Truncating.' if len(h_values) > n_points else 'Padding with last value.'}"
            )
            if len(h_values) > n_points:
                h_values = h_values[:n_points]
            else:
                pad = h_values[-1] if h_values else -1.0
                h_values.extend([pad] * (n_points - len(h_values)))

    # --- Build edge weight list ---
    # Format: for each node, for each neighbour index (0..deg), output the weight
    # if neighbour > node (undirected edge, written once per pair).
    # The reader expects one weight per undirected edge across all nodes.
    weight_entries = []
    for node in range(n_points):
        nbrs = neighbours(node)
        for i, nbr in enumerate(nbrs):
            if nbr > node:
                w = J_weights[node]  # could also vary per edge if desired
                weight_entries.append(w)

    # --- Write the file ---
    with open(output_path, "w") as f:
        f.write("$temp\n")
        f.write(f"{temp}\n\n")

        f.write("$dim_sizes\n")
        f.write(f"{dim}")
        for s in sizes:
            f.write(f" {s}")
        f.write("\n\n")

        f.write("$edge_weights_start\n")
        for w in weight_entries:
            f.write(f"{w}\n")
        f.write("\n")

        f.write("$mu\n")
        f.write(f"{mu}\n\n")

        f.write("$external_magnetic_field\n")
        for hv in h_values:
            f.write(f"{hv}\n")
        f.write("\n")

    print(f"Generated Ising model file: {output_path}")
    print(f"  Dimensions: {dim}D, sizes: {sizes} ({n_points} total spins)")
    print(f"  Couplings J:  {min(J_weights):.3f} to {max(J_weights):.3f} (all > 0 = ferromagnetic)")
    print(f"  Field h:      {min(h_values):.3f} to {max(h_values):.3f} (all < 0 = favours spin-down)")
    print(f"  Temperature:  {temp}")
    print(f"  Ground state: all zeros (all spins false / spin-down / σ = −1)")


def cli() -> None:
    parser = argparse.ArgumentParser(
        description="Generate an Ising model input file with ground state all-zeros."
    )
    parser.add_argument(
        "-d", "--dim", type=int, default=2,
        help="Number of dimensions (default: 2)"
    )
    parser.add_argument(
        "--sizes", type=int, nargs="+", default=[4, 4],
        help="Grid sizes per dimension (default: 4 4)"
    )
    parser.add_argument(
        "-t", "--temp", type=float, default=1.0,
        help="Temperature (default: 1.0)"
    )
    parser.add_argument(
        "-J", "--coupling", type=float, default=None,
        dest="J",
        help="Uniform ferromagnetic coupling weight (default: random 0.5–2.0 per edge)"
    )
    parser.add_argument(
        "-H", "--field", type=float, default=None,
        dest="h",
        help="Uniform external magnetic field per site (default: random −2.0 to −0.1 "
             "per site; negative favours spin-down)"
    )
    parser.add_argument(
        "-m", "--mu", type=float, default=1.0,
        help="Magnetic moment scaling factor (default: 1.0)"
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("ising_input.txt"),
        help="Output file path (default: ising_input.txt)"
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=None,
        help="Random seed for reproducible generation"
    )
    args = parser.parse_args()

    # Validate
    if args.dim < 1:
        parser.error("--dim must be >= 1")
    if len(args.sizes) != args.dim:
        parser.error(
            f"--sizes needs exactly {args.dim} values "
            f"(got {len(args.sizes)}: {args.sizes})"
        )
    if any(s < 1 for s in args.sizes):
        parser.error("All sizes must be >= 1")

    generate_ising_file(
        dim=args.dim,
        sizes=args.sizes,
        temp=args.temp,
        J=args.J,
        h=args.h,
        mu=args.mu,
        output_path=args.output,
        seed=args.seed,
    )


if __name__ == "__main__":
    cli()