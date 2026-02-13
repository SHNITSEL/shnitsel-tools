from .dim_red_result import DimRedResult
from .pca import PCA, pca, pca_direct, PCAResult
from .lda import LDA, lda, linear_discriminat_analysis, LDAResult
from .pls import PLS, pls, PLSResult

__all__ = [
    "PCA",
    "pca",
    "pca_direct",
    "LDA",
    "lda",
    "linear_discriminat_analysis",
    "PLS",
    "pls",
    "DimRedResult",
    "PCAResult",
    "LDAResult",
    "PLSResult",
]
