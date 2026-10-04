"""tools/rte -- TPRT-hallmark annotation of PEAR-TREE insertions (annotate_v2 plug-in).

Modules:
    library       reference RTE library + mappy indices (resources/rte_library layout)
    inputs        combine -> annotate sidecars (insertions.evidence.tsv.gz, insertions.reads.fa.gz)
    assembly      read layouts + covered-element consensus (indel-aware pileup)
    structure     element / 5' structure / tags
    transduction  3' transduction source lookup + novel-source rule
    pseudogene    exon-exon junction proof for processed pseudogenes
    hallmarks     TSD, EN motif, poly-A, slippage, fold-back
    score         transparent TPRT point system
    record        structured per-insertion record (annotate output columns)
    annotator     orchestration (RteAnnotator)
    calibrate     truth-vs-annotation feature separation / ROC
"""
