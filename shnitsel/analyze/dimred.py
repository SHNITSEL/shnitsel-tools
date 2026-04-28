from .dim_red_result import DimRedResult
from .pca import PCA, pca, pca_direct, PCAResult
from .lda import (
    LDA,
    lda,
    linear_discriminat_analysis,
    LDAResult,
    qda,
    QDA,
    quadratic_discriminat_analysis,
    QDAResult,
)
from .pls import PLS, pls, PLSResult

__all__ = [
    "PCA",
    "pca",
    "pca_direct",
    "LDA",
    "lda",
    "linear_discriminat_analysis",
    # "qda",
    # "QDA",
    # "quadratic_discriminat_analysis",
    # "QDAResult",
    "PLS",
    "pls",
    "DimRedResult",
    "PCAResult",
    "LDAResult",
    "PLSResult",
]
