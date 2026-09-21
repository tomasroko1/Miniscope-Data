"""
Módulo de análisis de co-activación y correlación neuronal.
"""

from .metrics import (
    compute_pairwise_correlation,
    compute_jaccard_cooccurrence,
    compute_cross_correlation_lags,
    compute_shuffled_significance,
    find_coactive_assemblies
)

__all__ = [
    'compute_pairwise_correlation',
    'compute_jaccard_cooccurrence',
    'compute_cross_correlation_lags',
    'compute_shuffled_significance',
    'find_coactive_assemblies'
]
