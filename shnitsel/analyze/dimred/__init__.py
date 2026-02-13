from .dim_red_result import DimRedResult
from .pca import PCA, pca, pca_direct
from .lda import LDA, lda, linear_discriminat_analysis
from .pls import PLS, pls

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
]
