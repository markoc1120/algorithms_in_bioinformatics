import pytest

from ..scripts.global_affine import pairwise_alignment_affine


@pytest.mark.parametrize(
    "seq1,seq2,expected_score,expected_alignments",
    [
        (
            "acgtgtcaacgt",
            "acgtcgtagcta",
            24,
            {("acgtgtcaacgt", "acgtcgtagcta")},
        ),
        (
            "aataat",
            "aagg",
            22,
            {("aataat", "aagg--"), ("aataat", "aa--gg"), ("aataat", "a--agg")},
        ),
        (
            "tccagaga",
            "tcgat",
            29,
            {("tccagaga", "tc---gat")},
        ),
        (
            "ggcctaaaggcgccggtctttcgtaccccaaaatctcggcattttaagataagtgagtgttgcgttacactagcgatctaccgcgtcttatacttaagcgtatgcccagatctgactaatcgtgcccccggattagacgggcttgatgggaaagaacagctcgtctgtttacgtataaacagaatcgcctgggttcgc",
            "gggctaaaggttagggtctttcacactaaagagtggtgcgtatcgtggctaatgtaccgcttctggtatcgtggcttacggccagacctacaagtactagacctgagaactaatcttgtcgagccttccattgagggtaatgggagagaacatcgagtcagaagttattcttgtttacgtagaatcgcctgggtccgc",
            395,
            {},
        ),
    ],
)
def test_pairwise_alignment_linear(seq1, seq2, expected_score, expected_alignments):
    S, alignments = pairwise_alignment_affine(seq1, seq2)
    assert S[-1, -1] == expected_score
    if not expected_alignments:
        return

    for alignment in alignments:
        assert alignment in expected_alignments
        expected_alignments.remove(alignment)
    assert not expected_alignments
