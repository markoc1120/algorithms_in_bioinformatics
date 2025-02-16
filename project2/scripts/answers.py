import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

from scripts.constants import GAPEXTEND, GAPOPEN
from scripts.global_affine import pairwise_alignment_affine
from scripts.global_linear import pairwise_alignment_linear
from scripts.helpers import format_seq, read_sequences


def parse_args():
    parser = argparse.ArgumentParser(description="Generate alignment analysis report")
    parser.add_argument("sequences", help="Path to input FASTA file")
    return parser.parse_args()


def format_questions(sequences: dict) -> str:
    seq1, seq2 = sequences["seq1"], sequences["seq2"]
    linear_cost, linear_alignments = pairwise_alignment_linear(seq1, seq2)
    affine_cost, affine_alignments = pairwise_alignment_affine(seq1, seq2)

    output = []
    output.append("### Question 1\n")
    output.append(f"Linear gap cost alignment (g(k)={GAPOPEN}*k):\n")
    output.append(f"Optimal score: {linear_cost[-1, -1]}\n")
    output.append("Optimal alignment:\n")
    align1, align2 = linear_alignments[0]
    output.append(
        f"```\nAligned sequence 1 ({len(align1)} bp):\n{format_seq(align1)}\n```"
    )
    output.append(
        f"```\nAligned sequence 2 ({len(align2)} bp):\n{format_seq(align2)}\n```"
    )

    output.append("### Question 2\n")
    output.append(f"Affine gap cost alignment (g(k)={GAPOPEN}+{GAPEXTEND}k):\n")
    output.append(f"Optimal score: {affine_cost[-1, -1]}\n")
    output.append("Optimal alignment:\n")
    align1, align2 = affine_alignments[0]
    output.append(
        f"```\nAligned sequence 1 ({len(align1)} bp):\n{format_seq(align1)}\n```"
    )
    output.append(
        f"```\nAligned sequence 2 ({len(align2)} bp):\n{format_seq(align2)}\n```"
    )

    output.append("### Question 3\n")
    output.append(f"Score matrix for linear gap cost (g(k)={GAPOPEN}k):\n")
    linear_matrix = create_score_matrix(sequences, "linear")
    output.append(f"```\n{linear_matrix}\n```\n")

    output.append("### Question 4\n")
    output.append(f"Score matrix for affine gap cost (g(k)={GAPOPEN}+{GAPEXTEND}k):\n")
    affine_matrix = create_score_matrix(sequences, "affine")
    output.append(f"```\n{affine_matrix}\n```\n")
    return "\n".join(output)


def create_score_matrix(sequences: dict, gap_model: str) -> pd.DataFrame:
    n = len(sequences)
    matrix = np.zeros((n, n))

    for (i, seq1), (j, seq2) in combinations(enumerate(sequences.values()), 2):
        if gap_model == "linear":
            score, _ = pairwise_alignment_linear(seq1, seq2)
        else:
            score, _ = pairwise_alignment_affine(seq1, seq2)
        matrix[i, j] = matrix[j, i] = score[-1, -1]

    return pd.DataFrame(
        matrix, index=sequences.keys(), columns=sequences.keys(), dtype=int
    )


def get_formatted_output(sequences_path: str) -> str:
    sequences = read_sequences(sequences_path)

    output = []
    output.append("## Sequence Alignment Analysiss\n")
    output.append(format_questions(sequences))

    return "".join(output)


if __name__ == "__main__":
    args = parse_args()
    formatted_output = get_formatted_output(args.sequences)
    print(formatted_output)
