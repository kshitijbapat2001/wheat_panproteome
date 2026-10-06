#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Single source of truth for table and figure numbering."""

MAIN_TABLES = [
    ("summary", None, "Summary of the pan-proteome with recruitment specificity"),
    ("locus_counts", "locus_counts_by_cultivar.csv", "Protein models and gene loci per cultivar"),
    ("cs_loci", "cs_locus_table.csv", "Reference gene loci underlying the proteotypes"),
    ("granularity", "granularity_partition.csv",
     "Prevalence partition at three clustering granularities"),
    ("rca_isoforms", "rca_isoform_summary.csv",
     "RCA proteotypes assigned to the nine established wheat RCA isoforms"),
    ("physicochemical", "physicochemical_features.csv",
     "Sequence-derived physicochemical descriptors"),
    ("candidates", "candidate_table.csv",
     "Prioritised candidate proteotypes for experimental follow-up"),
]

SUPP_TABLES = [
    ("diamond_hits_all", "diamond_hits_all.csv",
     "Complete DIAMOND homology search output"),
    ("cpn60_classes", "audit/cpn60_classes.csv",
     "Chaperonin subunit class assigned to every Cpn60 proteotype by annotation transfer"),
    ("proteotypes", "proteotypes.csv", "Non-redundant proteotype catalogue"),
    ("pav_matrix", "pav_matrix.csv", "Proteotype by cultivar presence-absence matrix"),
    ("novel_divergence", "audit/novel_divergence.csv",
     "Divergence of each reference-absent proteotype from its nearest reference sequence"),
    ("members", "members.csv",
     "Protein models resolved to gene loci, chromosomes and subgenomes"),
    ("rca_isoform_assignment", "audit/rca_isoform_assignment.csv",
     "Assignment of every RCA proteotype to the established wheat RCA isoforms"),
    ("proteotype_annotation", "proteotype_annotation.csv",
     "UniProtKB and InterPro annotation of every proteotype"),
    ("polymorphic_sites", "polymorphic_sites.csv",
     "Per-column polymorphism across the family alignments"),
    ("domain_classification", "domains/domain_classification.csv",
     "HMM domain profiling, Pfam confirmation and reciprocal best hits for every proteotype"),
    ("proteotype_recovery", "transcripts/proteotype_recovery.csv",
     "Recovery of each proteotype in an independent de novo annotation"),
    ("physicochemical_features", "physicochemical_features.csv",
     "Sequence-derived physicochemical descriptors for every proteotype"),
    ("proteotype_naming", "proteotype_naming.csv",
     "Composite identifiers and closest UniProt entry for every proteotype"),
    ("genome_wide_pav", "genome_wide_pav.csv",
     "Genome-wide proteotype-cluster presence-absence table"),
    ("pangenome_partition", "pangenome_partition.csv",
     "Core, soft-core, shell and cloud partition of the genome-wide clusters"),
    ("group_comparison", "group_comparison.csv",
     "Phenotype-stratified comparison of descriptors between cultivar groups "
     "(empty in this study: no phenotype labels were assigned)"),
]

SUPP_FIGURES = [
    ("RCA_PAV_heatmap", ["RCA_PAV_heatmap.png"], "Presence-absence matrix, RCA"),
    ("CPN60_PAV_heatmap", ["CPN60_PAV_heatmap.png"], "Presence-absence matrix, Cpn60"),
    ("RCA_alignment", ["RCA_alignment.pdf", "RCA_alignment_p1.png"],
     "Residue-level alignment of the RCA proteotypes"),
    ("CPN60_alignment", ["CPN60_alignment.pdf", "CPN60_alignment_p1.png"],
     "Residue-level alignment of the Cpn60 proteotypes by subunit class"),
    ("genomewide_prevalence_spectrum", ["genomewide_prevalence_spectrum.png"],
     "Prevalence spectrum of the genome-wide pan-proteome"),
    ("genomewide_per_cultivar", ["genomewide_per_cultivar.png"],
     "Per-cultivar composition of the genome-wide pan-proteome"),
]


def main_table_number(key: str) -> int:
    return [k for k, _, _ in MAIN_TABLES].index(key) + 1


def supp_table_number(key: str) -> int:
    return [k for k, _, _ in SUPP_TABLES].index(key) + 1


def supp_figure_number(key: str) -> int:
    return [k for k, _, _ in SUPP_FIGURES].index(key) + 1


SLABEL = {k: f"S{i}" for i, (k, _, _) in enumerate(SUPP_TABLES, 1)}
TLABEL = {k: str(i) for i, (k, _, _) in enumerate(MAIN_TABLES, 1)}
