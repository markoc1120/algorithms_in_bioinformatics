import os
import random

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .constants import (
    BASE_DECODE,
    CHUNK_SIZE,
    GAPEXTEND,
    GAPOPEN,
    SCORE_MATRIX_PATH,
    WIDTH,
)


def init_C(seq1_length: int, seq2_length: int) -> np.ndarray:
    C = np.zeros((seq1_length, seq2_length))
    C[0, :] = range(seq2_length)
    C[:, 0] = range(seq1_length)
    C[0, :] = C[0, :] * GAPOPEN
    C[:, 0] = C[:, 0] * GAPOPEN
    return C


def format_seq(seq_str: str) -> str:
    formatted_output = [
        char if idx % CHUNK_SIZE else char + "\n" for idx, char in enumerate(seq_str, 1)
    ]
    return "".join(formatted_output)


def format_output(seq1: str, seq2: str, costs: np.ndarray, results: list) -> str:
    separator = "-" * WIDTH + "\n"
    output = []
    output.append("```\n")
    output.append(separator)
    output.append("Alignment Plus Analysis\n")
    output.append(separator)
    output.append("Input Sequences:\n\n")
    output.append(f"Sequence 1 ({len(seq1)} bp):\n{format_seq(seq1)}\n")
    output.append(f"Sequence 2 ({len(seq2)} bp):\n{format_seq(seq2)}\n")
    output.append(separator)

    output.append("Alignment Statistics:\n\n")
    output.append(f"Gap open cost: {GAPOPEN}\n")
    output.append(f"Gap extend cost: {GAPEXTEND}\n")
    df = pd.read_csv(SCORE_MATRIX_PATH)
    df.index = BASE_DECODE.keys()
    output.append(f"Score matrix:\n{df}\n")
    output.append(f"Maximum alignment score: {int(costs[len(seq1), len(seq2)])}\n")
    output.append(separator)

    align1, align2 = results[random.randint(0, len(results) - 1)]
    output.append("Random Optimal Alignment:\n\n")
    output.append(f"Aligned sequence 1 ({len(seq1)} bp):\n{format_seq(align1)}\n")
    output.append(f"Aligned sequence 2 ({len(seq2)} bp):\n{format_seq(align2)}\n")
    output.append(separator)
    output.append("\n```")
    return "".join(output)


def read_sequence(path: str) -> str:
    try:
        return getattr(SeqIO.read(path, "fasta"), "seq", "")
    except FileNotFoundError:
        return path


def read_sequences(fasta_path: str) -> dict:
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}


def save_result(results: list[tuple[str, str]], output_path: str) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    records = []
    for idx, (seq1, seq2) in enumerate(results, 1):
        record1 = SeqRecord(
            Seq(seq1), id=f"alignment_{idx}_seq1", description="First sequence"
        )
        record2 = SeqRecord(
            Seq(seq2), id=f"alignment_{idx}_seq2", description="Second sequence"
        )
        records.extend([record1, record2])
    SeqIO.write(records, output_path, "fasta")
