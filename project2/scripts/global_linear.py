import numpy as np

from .constants import BASE_DECODE, GAPOPEN, SCORE_MATRIX
from .helpers import init_C


def pairwise_alignment_linear(
    seq1: str, seq2: str
) -> tuple[np.ndarray, list[tuple[str, str]]]:
    C = calc_cost_linear(seq1, seq2)
    alignments = backtrack_linear(C, seq1, seq2)
    return C, alignments


def calc_cost_linear(seq1: str, seq2: str) -> np.ndarray:
    # init C table based on sequences lengths
    seq1_length, seq2_length = len(seq1) + 1, len(seq2) + 1
    C = init_C(seq1_length, seq2_length)

    # fill C row by row
    for s1_idx in range(1, seq1_length):
        for s2_idx in range(1, seq2_length):
            match_cost = SCORE_MATRIX[
                BASE_DECODE[seq1[s1_idx - 1].upper()],
                BASE_DECODE[seq2[s2_idx - 1].upper()],
            ]
            # chosse minimum cost for the three cases
            C[s1_idx, s2_idx] = min(
                C[s1_idx - 1, s2_idx] + GAPOPEN,  # deletion
                C[s1_idx, s2_idx - 1] + GAPOPEN,  # insertion
                C[s1_idx - 1, s2_idx - 1] + match_cost,  # match/mismatch
            )
    return C


def backtrack_linear(C: np.ndarray, seq1: str, seq2: str) -> list[tuple[str, str]]:
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
            and C[s1_idx, s2_idx]
            == C[s1_idx - 1, s2_idx - 1]
            + SCORE_MATRIX[
                BASE_DECODE[seq1[s1_idx - 1].upper()],
                BASE_DECODE[seq2[s2_idx - 1].upper()],
            ]
        ):
            trace_alignments(
                s1_idx - 1,
                s2_idx - 1,
                seq1[s1_idx - 1] + align1,
                seq2[s2_idx - 1] + align2,
            )
        # check if we found a deletion -> add recursive call to stack
        if s1_idx > 0 and C[s1_idx - 1, s2_idx] + GAPOPEN == C[s1_idx, s2_idx]:
            trace_alignments(
                s1_idx - 1, s2_idx, seq1[s1_idx - 1] + align1, "-" + align2
            )
        # check if we found an insertion -> add recursive call to stack
        if s2_idx > 0 and C[s1_idx, s2_idx - 1] + GAPOPEN == C[s1_idx, s2_idx]:
            trace_alignments(
                s1_idx, s2_idx - 1, "-" + align1, seq2[s2_idx - 1] + align2
            )

    trace_alignments(s1_idx, s2_idx, align1, align2)
    return results
