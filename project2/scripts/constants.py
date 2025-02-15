import numpy as np
import pandas as pd


def read_gapcost(path):
    with open(path) as f:
        gapcost = f.read()
    return int(gapcost)


def read_score_matrix(path):
    df = pd.read_csv(path)
    return np.array(df.values)


BASE_DECODE = {"A": 0, "C": 1, "G": 2, "T": 3}
SCORE_MATRIX_PATH = "parameters/score_matrix.csv"
GAPOPEN_PATH = "parameters/gapopen"
GAPEXTEND_PATH = "parameters/gapextend"

SCORE_MATRIX = read_score_matrix(SCORE_MATRIX_PATH)
GAPOPEN = read_gapcost(GAPOPEN_PATH)  # 5
GAPEXTEND = read_gapcost(GAPEXTEND_PATH)  # 5
DEF_VAL = float("inf")

CHUNK_SIZE = 50
WIDTH = 80
