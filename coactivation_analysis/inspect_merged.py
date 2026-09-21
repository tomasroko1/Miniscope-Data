import argparse

import numpy as np
from scipy.stats import pearsonr


def main():
    parser = argparse.ArgumentParser(
        description="Compara dos primeras subsesiones de una matriz de co-activación."
    )
    parser.add_argument(
        "matrix_path",
        help="Ruta a corr_matrices_raw.npy generada por el análisis multi-subsesión.",
    )
    args = parser.parse_args()

    matrices = np.load(args.matrix_path)
    if matrices.ndim != 3 or matrices.shape[0] < 2:
        raise ValueError("La matriz debe tener al menos dos subsesiones.")

    triu_idx = np.triu_indices(matrices.shape[1], k=1)
    first = matrices[0][triu_idx]
    second = matrices[1][triu_idx]
    valid = np.isfinite(first) & np.isfinite(second)
    first = first[valid]
    second = second[valid]

    print(
        "Correlación global entre Sub 1 y Sub 2 "
        f"({len(first)} pares): {pearsonr(first, second)[0]:.3f}"
    )

    strong = first > 0.1
    print(f"Pares con r > 0.1 en Sub 1: {np.sum(strong)}")
    if np.any(strong):
        print(f"Media de r en Sub 1: {np.mean(first[strong]):.3f}")
        print(
            "Media de r en Sub 2 para esos mismos pares: "
            f"{np.mean(second[strong]):.3f}"
        )

    shuffled = np.random.default_rng().permutation(second)
    print(
        "Correlación nula entre matrices: "
        f"{pearsonr(first, shuffled)[0]:.3f}"
    )


if __name__ == "__main__":
    main()
