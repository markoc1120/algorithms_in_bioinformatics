import argparse
import sys

from scripts.global_affine import pairwise_alignment_affine
from scripts.global_linear import pairwise_alignment_linear
from scripts.helpers import format_output, save_result


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sequence alignment with linear and affine gap penalties."
    )
    parser.add_argument(
        "--seq1",
        type=str,
        required=True,
        help="First sequence to pairwise align, it can be either path of the sequence itself (required)",
    )
    parser.add_argument(
        "--seq2",
        type=str,
        required=True,
        help="Second sequence to pairwise align, it can be either path of the sequence itself (required)",
    )
    parser.add_argument(
        "--gap-model",
        choices=["linear", "affine"],
        default="linear",
        help="Gap penalty model to use, cost parameters can be modified in /parameters/... (default: linear)",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        help="Output path for the aligned file as a fasta, (e.g.: `--output-path results/alignments/XYZ.fasta`)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        if args.gap_model == "linear":
            costs, results = pairwise_alignment_linear(args.seq1, args.seq2)
        else:
            costs, results = pairwise_alignment_affine(args.seq1, args.seq2)

        if args.output_path:
            save_result(results, args.output_path)
        else:
            formatted_result = format_output(args.seq1, args.seq2, costs, results)
            print(formatted_result)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
