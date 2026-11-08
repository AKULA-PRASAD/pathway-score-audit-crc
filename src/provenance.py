"""Approved provenance statements for GSE39582.

These strings are APPROVED WORDING. They are written verbatim into every
analysis-bundle manifest. Do not paraphrase them: the qualifications they
carry are the point.
"""

PROCESSING_PROVENANCE = (
    "All 585 samples are annotated with the same deposited processing "
    "description, 'RMA normalisation (Bioconductor R affy package) + ComBat "
    "(R sva package)'. The deposit does not document the fitting population "
    "of either procedure - which arrays entered a given RMA or ComBat run is "
    "not stated. The 566 tumour samples and the 19 non-tumoural samples were "
    "submitted ~3 years apart (2012-07-23 and 2015-05-12) in separate "
    "accession blocks, and the ArrayExpress record carrying the ComBat batch "
    "assignments covers only the 566 tumour samples. It must therefore not be "
    "assumed that a single joint normalisation or batch-correction was fitted "
    "across every array. Expression values are used as deposited; this "
    "project performs no normalisation and no batch correction."
)

PARTITION_INDEPENDENCE = (
    "Discovery and CIT validation have distinct sample partitions but "
    "associated processing history. Batch composition differs sharply between "
    "them (batch01 33/0, batch02 14/0, batch09 1/14). Distinct partitions do "
    "NOT imply independence of processing-related information. CIT validation "
    "is internal validation only. Batch labels are preserved; no batch "
    "correction is applied and no batch covariate is automatically included."
)

TIME_UNIT_ASSUMPTION = (
    "months - PROVISIONAL working assumption, undocumented in any source. Six "
    "documentary sources were searched; none labels the rfs.delay column."
)

CENSORING_INTERPRETATION = (
    "death without recurrence = censoring - PROVISIONAL, inferred from a "
    "49/49 data pattern. NOT a published rule; the publication and its "
    "supplementary methods do not state this."
)

ENDPOINT_HIERARCHY = {
    "rfs_full": "primary",
    "rfs_60units": "conditional (depends on the provisional months assumption)",
    "os": "separate_preserved (not combined with RFS by Phase 3)",
}

# Phrases that must never appear in a manifest: they assert a fitting
# population the deposit does not document.
FORBIDDEN_PROVENANCE_PHRASES = (
    "jointly rma",
    "all 585 arrays",
    "across all 585",
    "jointly normalised across",
    "jointly normalized across",
)
