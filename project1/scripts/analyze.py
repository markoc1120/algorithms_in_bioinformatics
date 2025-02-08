import random
import sys

import numpy as np
from Bio import SeqIO

BASE_DECODE = {"A": 0, "C": 1, "G": 2, "T": 3}
SCORE_MATRIX = np.array([[10, 2, 5, 2], [2, 10, 2, 5], [5, 2, 10, 2], [2, 5, 2, 10]])
GAPCOST = -5


def read_sequences(path: str) -> str:
    try:
        return getattr(SeqIO.read(path, "fasta"), "seq", "")
    except FileNotFoundError:
        return path


def pairwise_alignment(seq1: str, seq2: str) -> None:
    # init dp table based on sequences lengths
    seq1_length, seq2_length = len(seq1) + 1, len(seq2) + 1
    dp = np.zeros((seq1_length, seq2_length))

    # init first row and column given the GAPCOST
    dp[0, :] = range(seq2_length)
    dp[:, 0] = range(seq1_length)
    dp[0, :] = dp[0, :] * GAPCOST
    dp[:, 0] = dp[:, 0] * GAPCOST

    # fill dp row by row
    for s1_idx in range(1, seq1_length):
        for s2_idx in range(1, seq2_length):
            match_cost = SCORE_MATRIX[
                BASE_DECODE[seq1[s1_idx - 1]], BASE_DECODE[seq2[s2_idx - 1]]
            ]
            # chosse minimum cost for the three cases
            dp[s1_idx, s2_idx] = max(
                dp[s1_idx - 1, s2_idx] + GAPCOST,  # deletion
                dp[s1_idx, s2_idx - 1] + GAPCOST,  # insertion
                dp[s1_idx - 1, s2_idx - 1] + match_cost,  # match/mismatch
            )
    alignments = backtrack(dp, seq1, seq2)
    return dp, alignments


def backtrack(dp: np.ndarray, seq1: str, seq2: str) -> list[tuple[str, str]]:
    s1_idx, s2_idx = len(seq1), len(seq2)
    align1 = align2 = ""
    results = []

    def trace_alignments(s1_idx, s2_idx, align1, align2):
        # check if we have depleted our indexes, save the two alignments
        if s1_idx == s2_idx == 0:
            results.append((align1, align2))
            return

        # check if we found a match/mismatch -> add recursive call to stack
        if (
            s1_idx > 0
            and s2_idx > 0
            and dp[s1_idx, s2_idx]
            == dp[s1_idx - 1, s2_idx - 1]
            + SCORE_MATRIX[BASE_DECODE[seq1[s1_idx - 1]], BASE_DECODE[seq2[s2_idx - 1]]]
        ):
            trace_alignments(
                s1_idx - 1,
                s2_idx - 1,
                seq1[s1_idx - 1] + align1,
                seq2[s2_idx - 1] + align2,
            )
        # check if we found a deletion -> add recursive call to stack
        if s1_idx > 0 and dp[s1_idx - 1, s2_idx] + GAPCOST == dp[s1_idx, s2_idx]:
            trace_alignments(
                s1_idx - 1, s2_idx, seq1[s1_idx - 1] + align1, "-" + align2
            )
        # check if we found an insertion -> add recursive call to stack
        if s2_idx > 0 and dp[s1_idx, s2_idx - 1] + GAPCOST == dp[s1_idx, s2_idx]:
            trace_alignments(
                s1_idx, s2_idx - 1, "-" + align1, seq2[s2_idx - 1] + align2
            )

    trace_alignments(s1_idx, s2_idx, align1, align2)
    return results


def count_alignments(n: int, m: int, memo: dict[tuple[int, int] : int]) -> int:
    key = (n, m)
    if key in memo:
        return memo[key]

    # base cases
    if n <= 0 or m <= 0:
        return 1

    count = (
        count_alignments(n - 1, m - 1, memo)  # previous diagonal
        + count_alignments(n - 1, m, memo)  # previous vertical
        + count_alignments(n, m - 1, memo)  # previous horizontal
    )
    memo[key] = count
    return count


def format_output(
    seq1: str, seq2: str, dp: np.ndarray, results: list, upper_bound: int
) -> str:
    width, chunk_size = 80, 50
    separator = "-" * width + "\n"

    output = []
    output.append("```\n")
    output.append(separator)
    output.append("Sequence Alignment Analysis\n")
    output.append(separator)

    output.append("Input Sequences:\n")
    output.append(
        f"Sequence 1 ({len(seq1)} bp): {seq1[:chunk_size]}{'...' if len(seq1) > chunk_size else ''}\n"
    )
    output.append(
        f"Sequence 2 ({len(seq2)} bp): {seq2[:chunk_size]}{'...' if len(seq2) > chunk_size else ''}\n"
    )
    output.append(separator)

    output.append("Alignment Statistics:\n")
    output.append(f"Maximum alignment score: {int(dp[len(seq1), len(seq2)])}\n")
    if upper_bound > 1e10:
        output.append(f"Total possible alignments: {upper_bound:.2e}\n")
    else:
        output.append(f"Total possible alignments: {upper_bound:,}\n")
    output.append(f"Number of optimal alignments: {len(results):,}\n")
    output.append(separator)

    align1, align2 = results[random.randint(0, len(results) - 1)]
    output.append("Random Optimal Alignment:\n")

    output.append(
        f"\nSeq1: {align1[:chunk_size]}{'...' if len(align1) > chunk_size else ''}\n"
    )
    output.append(
        f"Seq2: {align2[:chunk_size]}{'...' if len(align2) > chunk_size else ''}\n"
    )
    output.append(separator)
    output.append("```\n")
    output.append("\pagebreak")
    return "".join(output)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python main.py /path/to/seq1 /path/to/seq1")
        sys.exit(1)

    seq1_path, seq2_path = sys.argv[1], sys.argv[2]
    try:
        seq1, seq2 = read_sequences(seq1_path), read_sequences(seq2_path)
    except Exception:
        print("There was a problem while reading fasta files")
        sys.exit(1)

    dp, results = pairwise_alignment(seq1, seq2)
    upper_bound = count_alignments(len(seq1), len(seq2), {})
    formatted_output = format_output(seq1, seq2, dp, results, upper_bound)
    print(formatted_output)
