"""pathway-audit: audit what a gene-set score measures in bulk tumour expression data.

Core idea: before interpreting a pathway score as biology, check (1) that the expression data pass positive-control
co-expression QC, (2) which genes drive the score and whether it is internally coherent, (3) how much of it is tissue
composition, (4) whether genotype associations survive adjustment for composition and covariates, and optionally
(5) whether it adds optimism-corrected discrimination to a clinical survival model.
"""
__version__ = "0.1.0"

from .audit import audit  # noqa: F401
from .qc import coexpression_qc  # noqa: F401
from .scoring import score_gene_set, zscore_rows  # noqa: F401
