import pytest

from ..scripts.global_linear import pairwise_alignment_linear


@pytest.mark.parametrize(
    "seq1,seq2,expected_score,expected_alignments",
    [
        (
            "acgtgtcaacgt",
            "acgtcgtagcta",
            22,
            {("acgt-gtcaacgt-", "acgtcgt-agc-ta"), ("acgt-gtcaacgt", "acgtcgt-agcta")},
        ),
        ("aataat", "aagg", 14, {("aataat", "aa-gg-")}),
        (
            "tccagaga",
            "tcgat",
            20,
            {
                ("tccagaga", "tc--gat-"),
                ("tccagaga", "t-c-gat-"),
                ("tccagaga", "tc--ga-t"),
                ("tccagaga", "t-c-ga-t"),
            },
        ),
        (
            "ggcctaaaggcgccggtctttcgtaccccaaaatctcggcattttaagataagtgagtgttgcgttacactagcgatctaccgcgtcttatacttaagcgtatgcccagatctgactaatcgtgcccccggattagacgggcttgatgggaaagaacagctcgtctgtttacgtataaacagaatcgcctgggttcgc",
            "gggctaaaggttagggtctttcacactaaagagtggtgcgtatcgtggctaatgtaccgcttctggtatcgtggcttacggccagacctacaagtactagacctgagaactaatcttgtcgagccttccattgagggtaatgggagagaacatcgagtcagaagttattcttgtttacgtagaatcgcctgggtccgc",
            325,
            {},
        ),
    ],
)
def test_pairwise_alignment_linear(seq1, seq2, expected_score, expected_alignments):
    C, alignments = pairwise_alignment_linear(seq1, seq2)
    assert C[-1, -1] == expected_score
    if not expected_alignments:
        return

    for alignment in alignments:
        assert alignment in expected_alignments
        expected_alignments.remove(alignment)
    assert not expected_alignments
