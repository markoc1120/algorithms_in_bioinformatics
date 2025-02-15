import numpy as np

from .constants import BASE_DECODE, GAPEXTEND, GAPOPEN, SCORE_MATRIX


def pairwise_alignment_affine(
    seq1: str, seq2: str
) -> tuple[np.ndarray, list[tuple[str, str]]]:
    S = calc_cost_affine(seq1, seq2)
    alignments = backtrack_affine(S, seq1, seq2)
    return S, alignments


def calc_cost_affine(seq1: str, seq2: str) -> np.ndarray:
    seq1_length, seq2_length = len(seq1) + 1, len(seq2) + 1

    S = np.full((seq1_length, seq2_length), float("inf"))
    D = np.full((seq1_length, seq2_length), float("inf"))
    I = np.full((seq1_length, seq2_length), float("inf"))
    S[0, 0] = 0

    # fill S, D, I row by row
    for i in range(seq1_length):
        for j in range(seq2_length):
            if i == j == 0:
                continue

            if i > 0 and j > 0:
                match_score = SCORE_MATRIX[
                    BASE_DECODE[seq1[i - 1].upper()], BASE_DECODE[seq2[j - 1].upper()]
                ]
                s_match = S[i - 1, j - 1] + match_score

            # calc D(i,j)
            if i > 0:
                d_scores = []
                if j >= 0:
                    d_scores.extend(
                        [S[i - 1, j] + GAPOPEN + GAPEXTEND, D[i - 1, j] + GAPEXTEND]
                    )
                if d_scores:
                    D[i, j] = min(d_scores)

            # calc I(i,j)
            if j > 0:
                i_scores = []
                if i >= 0:
                    i_scores.extend(
                        [S[i, j - 1] + GAPOPEN + GAPEXTEND, I[i, j - 1] + GAPEXTEND]
                    )
                if i_scores:
                    I[i, j] = min(i_scores)

            # calc S(i,j)
            s_scores = []
            if i > 0 and j > 0:
                s_scores.append(s_match)
            if i > 0 and j >= 0:
                s_scores.append(D[i, j])
            if i >= 0 and j > 0:
                s_scores.append(I[i, j])
            if s_scores:
                S[i, j] = min(s_scores)
    return S


def backtrack_affine(S: np.ndarray, seq1: str, seq2: str) -> list[tuple[str, str]]:
    s1_idx, s2_idx = len(seq1), len(seq2)
    align1 = align2 = ""
    results = []

    def trace_alignments(s1_idx, s2_idx, align1, align2):
        # check if we have depleted our indexes, save the two alignments
        if s1_idx == s2_idx == 0:
            results.append((align1, align2))
            return

        curr_score = S[s1_idx, s2_idx]
        # check if we found a match/mismatch -> add recursive call to stack
        if (
            s1_idx > 0
            and s2_idx > 0
            and curr_score
            == S[s1_idx - 1, s2_idx - 1]
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

        # try deletion(s) and insertion(s)
        k = 1
        while s1_idx >= k or s2_idx >= k:
            # check if we found deletion(s) -> add recursive call to stack
            if s1_idx >= k and curr_score == S[s1_idx - k, s2_idx] + GAPOPEN + (
                k * GAPEXTEND
            ):
                trace_alignments(
                    s1_idx - k,
                    s2_idx,
                    seq1[s1_idx - k : s1_idx] + align1,
                    "-" * k + align2,
                )

            # check if we found insertion(s) -> add recursive call to stack
            if s2_idx >= k and curr_score == S[s1_idx, s2_idx - k] + GAPOPEN + (
                k * GAPEXTEND
            ):
                trace_alignments(
                    s1_idx,
                    s2_idx - k,
                    "-" * k + align1,
                    seq2[s2_idx - k : s2_idx] + align2,
                )

            k += 1

    trace_alignments(s1_idx, s2_idx, align1, align2)
    return results
