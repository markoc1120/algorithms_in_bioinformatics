import argparse
import functools
import random
import timeit
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from scripts.global_affine import calc_cost_affine
from scripts.global_linear import calc_cost_linear

plt.rcParams["figure.figsize"] = [12, 8]


def parse_args():
    parser = argparse.ArgumentParser(description="Generate benchmark figure")
    parser.add_argument("--output-path", help="Output path for the figure")
    return parser.parse_args()


def generate_random_sequence(length: int) -> str:
    return "".join(random.choice("ATGC") for _ in range(length))


def run_benchmark(seq_lengths: np.ndarray) -> Tuple[list, list]:
    linear_times, affine_times = [], []

    for n in seq_lengths:
        seq1 = generate_random_sequence(n)
        seq2 = generate_random_sequence(n)

        linear_time = timeit.timeit(
            functools.partial(calc_cost_linear, seq1, seq2),
            number=5,
        )
        linear_times.append(linear_time)

        affine_time = timeit.timeit(
            functools.partial(calc_cost_affine, seq1, seq2),
            number=1,
        )
        affine_times.append(affine_time)
        print(f"Completed length {n}")

    return linear_times, affine_times


def plot_results(
    seq_lengths: np.ndarray, linear_times: list, affine_times: list, output_path: str
) -> None:
    plt.figure(figsize=(12, 8))

    plt.plot(seq_lengths, linear_times, "bo-", label="Linear Gap", alpha=0.5)
    plt.plot(seq_lengths, affine_times, "ro-", label="Affine Gap", alpha=0.5)

    degree = 2  # expected O(n²) complexity
    linear_coeffs = np.polyfit(seq_lengths, linear_times, degree)
    affine_coeffs = np.polyfit(seq_lengths, affine_times, degree)

    linear_fit = np.poly1d(linear_coeffs)
    affine_fit = np.poly1d(affine_coeffs)

    plt.plot(seq_lengths, [linear_fit(n) for n in seq_lengths], "b--", alpha=0.3)
    plt.plot(seq_lengths, [affine_fit(n) for n in seq_lengths], "r--", alpha=0.3)

    plt.xlabel("Sequence Length")
    plt.ylabel("Time (seconds)")
    plt.title("Performance Comparison: Linear vs Affine Gap Penalties")
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.savefig(output_path)


if __name__ == "__main__":
    seq_lengths = np.linspace(0, 1000, 8, dtype=int)
    linear_times, affine_times = run_benchmark(seq_lengths)
    args = parse_args()
    print(args.output_path)
    plot_results(seq_lengths, linear_times, affine_times, args.output_path)
