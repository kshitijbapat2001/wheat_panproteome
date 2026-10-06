#!/usr/bin/env python3
"""Format supplementary tables to publication standard."""
import csv, glob, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numbering

CITE = ("Bapat et al., 2026. A pan-proteome catalogue of bread wheat: Rubisco activase "
        "and chaperonin-60 as resources for dissecting the genetic basis of "
        "photosynthetic thermotolerance")

CULT = ("One column per cultivar. 1 indicates the cultivar contributes a protein exactly "
        "identical to this proteotype; 0 indicates it does not. A zero denotes sequence "
        "difference, not absence of the gene: RCA gene content is invariant at nine loci "
        "in all ten cultivars.")
CULT_GW = ("One column per cultivar. 1 indicates the cultivar contributes at least one "
           "protein to this genome-wide cluster; 0 indicates it does not.")
CULTS = ["chinese_spring","arinalrfor","jagger","julius","lancer","landmark","mace",
         "norin61","stanley","mattis"]

COMMON = {
 "proteotype_id":"Stable within-study identifier for one distinct amino-acid sequence, formatted <FAMILY>_PT<nnn>.",
 "family":"Protein family: RCA (Rubisco activase) or CPN60 (chaperonin-60).",
 "n_cultivars":"Number of the ten cultivars carrying a protein exactly identical to this proteotype.",
 "cultivars":"Semicolon-separated list of those cultivars.",
 "novel_vs_reference":"TRUE if no exactly identical sequence occurs in the Chinese Spring reference proteome. Because proteotypes are defined by exact identity, a single substitution sets this flag; the magnitude of divergence is quantified in Supplemental Table {novel_divergence}.",
 "source_ids":"Semicolon-separated list of every contributing protein model (transcript identifier) that collapsed into this proteotype.",
 "length":"Length in amino-acid residues.",
 "uniprot_id":"Accession of the closest UniProtKB entry, taken as the best DIAMOND hit across the proteotype's member models.",
 "uniprot_entry":"UniProtKB entry name corresponding to that accession.",
 "representative_cultivar":"Cultivar used in the composite identifier. Chinese Spring where the proteotype occurs in the reference; otherwise the member with the best e-value to the reference query, ties broken by lexicographic order of cultivar then source identifier so the choice is reproducible across builds.",
 "representative_source_id":"Protein model of that representative cultivar.",
 "composite_id":"Publication identifier in the form <proteotype>__<representative cultivar>__<UniProt accession>. Used in Table {T_candidates}, in figure row labels and as the FASTA header of the deposited catalogue.",
}

TABLES = {
 1:("diamond_hits_all","Complete DIAMOND homology search output",{
   "cultivar":"Cultivar whose annotated proteome was searched.",
   "qseqid":"Query identifier: the curated UniProt reference sequence used to search.",
   "sseqid":"Subject identifier: the cultivar protein model returned as a hit.",
   "pident":"Percentage amino-acid identity over the aligned region.",
   "length":"Alignment length in residues.",
   "mismatch":"Number of mismatched positions in the alignment.",
   "gapopen":"Number of gap openings in the alignment.",
   "qstart":"Start position of the alignment in the query sequence.",
   "qend":"End position of the alignment in the query sequence.",
   "sstart":"Start position of the alignment in the subject sequence.",
   "send":"End position of the alignment in the subject sequence.",
   "evalue":"Expectation value of the alignment.",
   "bitscore":"Bit score of the alignment.",
   "qcovhsp":"Percentage of the query sequence covered by the alignment.",
   "scovhsp":"Percentage of the subject sequence covered by the alignment."}),
 2:("cpn60_classes","Chaperonin subunit class assigned to every Cpn60 proteotype by annotation transfer",{
   "proteotype_id":COMMON["proteotype_id"],
   "cpn60_class":"Subunit class transferred from the annotation of the best-matching curated reference query: plastid_Cpn60alpha, plastid_Cpn60beta or mitochondrial_groupI. Partitions the 141 Cpn60 proteotypes as 35, 75 and 31 respectively.",
   "best_reference":"Curated UniProt reference query from which the class was transferred.",
   "identity_to_best_reference":"Percentage amino-acid identity to that reference query. Reported because class transfer by identity is weaker evidence than domain-level assignment; Supplemental Table {domain_classification} reaches the same partition independently."}),
 3:("proteotype_annotation","UniProtKB and InterPro annotation of every proteotype",{
   **{k:COMMON[k] for k in ["proteotype_id","family","uniprot_id","uniprot_entry","representative_cultivar","representative_source_id","n_cultivars","novel_vs_reference","composite_id","cultivars","source_ids"]},
   "protein_name":"Recommended protein name of the closest UniProtKB entry.",
   "gene_names":"Gene name or names recorded for that entry, semicolon-separated.",
   "subcellular_location":"Annotated subcellular location or locations, semicolon-separated. Chloroplast stroma for plastid members; supports the compartment assignment of the group-2 RCA-like clade.",
   "pfam":"Pfam domain accessions with entry names, semicolon-separated. PF00118 identifies group I chaperonins; PF00004 and PF21228 the AAA core and the RCA-diagnostic helical domain.",
   "interpro":"InterPro entry accessions with names, semicolon-separated. IPR044960 denotes Rubisco activase-like entries.",
   "panther":"PANTHER family accessions, semicolon-separated.",
   "sequence_length":"Length of the UniProtKB entry sequence in residues, for comparison with the proteotype length.",
   "reviewed":"TRUE if the UniProtKB entry is Swiss-Prot reviewed rather than TrEMBL.",
   "interpro_full":"Complete list of InterPro entries matched to the accession by the InterPro API, semicolon-separated."}),
 4:("novel_divergence","Divergence of each reference-absent proteotype from its nearest reference sequence",{
   "proteotype_id":"Reference-absent proteotype identifier.",
   "family":COMMON["family"],
   "length":COMMON["length"],
   "n_cultivars":COMMON["n_cultivars"],
   "nearest_reference_proteotype":"Reference-present proteotype of the same family to which it is most similar by global pairwise alignment.",
   "identity_pct":"Percentage identity to that nearest reference proteotype. Needleman-Wunsch global alignment, BLOSUM62, gap open -11, gap extend -1.",
   "residue_differences":"Number of differing residue positions in that alignment. Reported so that a form differing by a single residue is not read as equivalent to a substantially divergent protein."}),
 5:("physicochemical_features","Sequence-derived physicochemical descriptors for every proteotype",{
   "length":COMMON["length"],
   "length_std_aa":"Length after removal of non-standard residue codes; the length used for all descriptor calculations.",
   "mw":"Molecular weight in daltons, computed from amino-acid composition.",
   "pI":"Theoretical isoelectric point.",
   "aromaticity":"Combined frequency of phenylalanine, tryptophan and tyrosine.",
   "gravy":"Grand average of hydropathy (Kyte-Doolittle); higher values indicate greater overall hydrophobicity.",
   "aliphatic_index":"Relative volume occupied by aliphatic side chains; higher values are associated with greater thermostability in globular proteins.",
   "frac_acidic_DE":"Fraction of residues that are aspartate or glutamate.",
   "frac_basic_KR":"Fraction of residues that are lysine or arginine.",
   "frac_charged_DEKRH":"Fraction of residues that are charged (D, E, K, R, H).",
   "net_charge_frac":"Net charge (basic minus acidic residues) as a fraction of sequence length.",
   "frac_hydrophobic":"Fraction of residues in the hydrophobic set (A, V, I, L, M, F, W, C).",
   "frac_proline":"Fraction of residues that are proline; relevant to backbone rigidity.",
   "frac_glycine":"Fraction of residues that are glycine; relevant to backbone flexibility.",
   "frac_arg_of_RK":"Arginine as a proportion of all arginine plus lysine residues; arginine enrichment at the expense of lysine is a recognised thermostability signature.",
   "instability_index":"Guruprasad instability index; values above 40 predict a protein unstable in vitro.",
   "frac_helix":"Predicted fraction of residues in helical secondary structure.",
   "frac_sheet":"Predicted fraction of residues in beta-sheet secondary structure.",
   "proteotype_id":COMMON["proteotype_id"],
   "family":COMMON["family"],
   "novel_vs_reference":COMMON["novel_vs_reference"]}),
 6:("polymorphic_sites","Per-column polymorphism across the family alignments",{
   "family":"Protein family whose alignment the column belongs to: RCA or CPN60.",
   "aln_column":"One-based column index in that family's MAFFT alignment of all proteotypes. Alignment coordinates, not positions in any single protein.",
   "n_variants":"Number of distinct amino-acid residues observed in the column, excluding gaps.",
   "residues":"Observed residues with their counts, semicolon-separated.",
   "entropy_bits":"Shannon entropy of the residue distribution in the column, in bits. Zero indicates an invariant column. Plotted in Figures 5 and 6.",
   "gap_fraction":"Fraction of sequences carrying a gap at this column. High values indicate columns spanned by indels or partial gene models rather than by substitution."}),
 7:("genome_wide_pav","Genome-wide proteotype-cluster presence-absence table",{
   "cluster_id":"Identifier of a genome-wide proteotype cluster. Clusters anchored on a Chinese Spring protein carry that protein's identifier; clusters formed from proteins with no acceptable reference match are numbered NOVEL_nnnnnn.",
   **{c:CULT_GW for c in CULTS},
   "novel_vs_reference":"TRUE for clusters formed among reference-absent proteins, FALSE for clusters anchored on a Chinese Spring protein."}),
 8:("pangenome_partition","Core, soft-core, shell and cloud partition of the genome-wide clusters",{
   "category":"Prevalence class of a genome-wide proteotype cluster: core (present in all ten cultivars), soft-core (present in at least 90 percent), shell (present in more than one but fewer than 90 percent) or cloud (present in a single cultivar).",
   "n_clusters":"Number of clusters in that class. Totals 131,282 clusters, of which 71,974 (54.8 percent) are core. Computed on protein sequences with reference-anchored clustering and therefore not directly equivalent to gene-level pan-genome estimates."}),
 9:("group_comparison","Phenotype-stratified comparison of descriptors between cultivar groups",{
   "family":COMMON["family"],
   "feature":"Physicochemical descriptor tested, as defined in Supplemental Table {physicochemical_features}.",
   "n_resistant":"Number of proteotypes contributed by cultivars assigned to the heat-resistant group.",
   "n_sensitive":"Number of proteotypes contributed by cultivars assigned to the heat-sensitive group.",
   "median_resistant":"Median value of the descriptor in the resistant group.",
   "median_sensitive":"Median value of the descriptor in the sensitive group.",
   "U":"Mann-Whitney U statistic for the two-sided test between groups.",
   "p_value":"Uncorrected two-sided p-value.",
   "q_value_BH":"p-value adjusted for multiple testing across descriptors by the Benjamini-Hochberg procedure."}),
 10:("proteotypes","Non-redundant proteotype catalogue",{
   "proteotype_id":COMMON["proteotype_id"],"family":COMMON["family"],"length":COMMON["length"],
   "n_cultivars":COMMON["n_cultivars"],"cultivars":COMMON["cultivars"],
   "present_in_reference":"TRUE if an exactly identical sequence occurs in the Chinese Spring reference proteome.",
   "novel_vs_reference":COMMON["novel_vs_reference"],"source_ids":COMMON["source_ids"]}),
 11:("proteotype_naming","Composite identifiers and closest UniProt entry for every proteotype",{
   k:COMMON[k] for k in ["proteotype_id","family","uniprot_id","uniprot_entry","representative_cultivar","representative_source_id","n_cultivars","novel_vs_reference","composite_id","cultivars","source_ids"]}),
 12:("members","Protein models resolved to gene loci, chromosomes and subgenomes",{
   "cultivar":"Cultivar whose annotated proteome contributed this protein model.",
   "family":COMMON["family"],
   "source_id":"Ensembl protein identifier of the individual protein model. Ensembl proteomes carry one record per annotated transcript, so several rows may share a gene.",
   "gene":"Parent gene locus of the protein model, obtained by stripping the transcript suffix. This is the level at which gene content is counted.",
   "chromosome":"Chromosome carrying the gene locus; 'un' denotes an unanchored scaffold.",
   "subgenome":"Subgenome assignment derived from the chromosome name: A, B, D, or 'un'.",
   "homoeologous_group":"Homoeologous group (1-7) derived from the chromosome name, or 'un'. RCA loci fall into groups 2 and 4.",
   "pident_to_ref":"Percentage amino-acid identity of this protein model to its best-matching curated reference query.",
   "evalue_to_ref":"E-value of that best DIAMOND hit.",
   "length":COMMON["length"],
   "proteotype_id":"Proteotype into which this model collapsed."}),
 13:("pav_matrix","Proteotype by cultivar presence-absence matrix",{
   "proteotype_id":COMMON["proteotype_id"],"family":COMMON["family"],
   **{c:CULT for c in CULTS}}),
 14:("domain_classification","HMM domain profiling, Pfam confirmation and reciprocal best hits for every proteotype",{
   "proteotype_id":COMMON["proteotype_id"],
   "hmm_class":"Class assigned by the best-scoring class profile: plastid_Cpn60alpha, plastid_Cpn60beta, mitochondrial_CPN60, rubisco_activase or cbbX_negative_control. Profiles were built with HMMER from reviewed Swiss-Prot sequences.",
   "hmm_bitscore":"Bit score against the best-scoring class profile.",
   "hmm_evalue":"E-value of that best profile match.",
   "runner_up":"Second-best-scoring class profile.",
   "bitscore_margin":"Bit-score difference between the best and second-best class. A margin above 20 bits is required for a high-confidence call.",
   "rbh_partner":"Curated reference sequence returned as the best DIAMOND hit for this proteotype.",
   "rbh_class":"Class of that reference sequence, determined independently of the HMM assignment.",
   "reciprocal":"TRUE if the reference sequence also returns this proteotype as its own best hit. A strict one-to-one test, and therefore conservative where many paralogues compete for few reference sequences.",
   "has_PF00118":"TRUE if Pfam PF00118 (Cpn60_TCP1) is detected; diagnostic of group I chaperonins. Detected in all 141 Cpn60 proteotypes.",
   "has_PF00004":"TRUE if Pfam PF00004 (AAA ATPase core) is detected.",
   "has_PF21228":"TRUE if Pfam PF21228 (RuBisCO activase AAA helical domain) is detected; RCA-diagnostic. Detected in 45 of 46 RCA proteotypes.",
   "confidence":"Overall call. High confidence requires the profile assignment and the reciprocal best hit to agree on class and the bit-score margin to exceed 20 bits; all other outcomes are flagged rather than hidden. A CbbX profile was included as an out-of-family negative control and was the best-scoring class for no proteotype, demonstrating that the classification is capable of returning a negative result."}),
 15:("proteotype_recovery","Recovery of each proteotype in an independent de novo annotation",{
   "proteotype_id":COMMON["proteotype_id"],"family":COMMON["family"],
   "n_loci":"Number of distinct gene loci contributing protein models to this proteotype.",
   "n_recovered":"Number of those loci whose gene identifier also appears in the independent de novo annotation of White et al. (2025).",
   "genes":"Semicolon-separated list of the contributing gene identifiers.",
   "any_locus_recovered":"TRUE if at least one contributing locus is recovered by the independent annotation.",
   "fraction_recovered":"n_recovered divided by n_loci.",
   "novel":"TRUE if no exactly identical sequence occurs in Chinese Spring. Recovery is reported separately for reference-present and reference-absent proteotypes, since it is the status of the latter that the comparison is intended to test.",
   "n_cultivars":COMMON["n_cultivars"]}),
 "rca":("rca_isoform_assignment","Assignment of every RCA proteotype to the established wheat RCA isoforms",{
   "proteotype_id":COMMON["proteotype_id"],
   "composite_id":COMMON["composite_id"],
   "established_isoform":"Established wheat RCA isoform, named isoform-homoeologue: RCA1\u03b2, RCA2\u03b1 or RCA2\u03b2 encoded by the 4A, 4B or 4D homoeologue (e.g. RCA2\u03b1-4D), or RCA-like-2A, -2B or -2D for the divergent group 2 lineage, which is not one of the nine established isoforms. Where an identical sequence is encoded by two homoeologues both are given (e.g. RCA2\u03b1-4A/4D).",
   "gene":"Gene of the tandem Rca pair: TaRca1, which encodes only RCA1\u03b2, or TaRca2, which encodes RCA2\u03b1 and RCA2\u03b2 by alternative splicing; RCA-like (group 2) for the divergent lineage.",
   "isoform":"Isoform class: RCA1\u03b2; RCA2\u03b1, which carries the redox-regulated C-terminal extension; RCA2\u03b2; or RCA-like.",
   "chromosomes":"Chromosome or chromosomes of the parent gene loci.",
   "gene_loci":"Parent gene loci of every member protein model, across all cultivars, semicolon-separated.",
   "length":COMMON["length"],
   "n_cultivars":COMMON["n_cultivars"],
   "cultivars":COMMON["cultivars"],
   "novel_vs_reference":COMMON["novel_vs_reference"],
   "model_status":"full-length if the proteotype covers at least 97 percent of the canonical form of its isoform with no insertion or deletion longer than four residues; otherwise 'variant model', followed by the missing or additional segments in the coordinates of the canonical precursor (an internal deletion is typical of an alternative splice model).",
   "coverage_of_canonical_pct":"Percentage of the canonical form of the isoform (its most prevalent proteotype) covered by the alignment.",
   "canonical_proteotype":"Most prevalent full-length proteotype of the same isoform and homoeologue, used as the reference for substitutions; for example RCA_PT001 for RCA2\u03b1-4A.",
   "identity_to_canonical_pct":"Percentage identity to that canonical proteotype over the aligned residues.",
   "substitutions_vs_canonical":"Number of amino-acid substitutions relative to the canonical proteotype. Insertions and deletions are described under model_status, not counted here.",
   "substitution_list":"The substitutions, each written as canonical residue, position in the canonical precursor and proteotype residue (e.g. L260F), semicolon-separated.",
   "switch_site_M159I":"Residue at the thermal and regulatory switch of Degen et al. (2020): Met-159 in their numbering of RCA2\u03b2, residue 206 of the precursor. The M159I substitution raises the temperature optimum of RCA2\u03b2 by 5 degrees C. Met in every RCA2 proteotype and Ile in every RCA1\u03b2 proteotype in this panel.",
   "site_K428":"Residue at Lys-428 of RCA2\u03b1, the first residue of the C-terminal extension; its replacement by Arg greatly reduces ADP inhibition (Scafaro et al. 2019b). Reported for RCA2\u03b1 proteotypes only.",
   "cte_cysteine_pair":"Residues at the two cysteines of the C-terminal extension (Cys-441 and Cys-460 of the RCA2\u03b1 precursor), the redox-regulated pair reduced by thioredoxin-f (Zhang & Portis 1999). Reported for RCA2\u03b1 proteotypes only; a dash marks a position absent from the model.",
   "cte_residues_aligned":"Number of residues of the RCA2\u03b1 C-terminal extension (residues 428-464 of the canonical RCA2\u03b1) found in the proteotype by local alignment at 60 percent identity or more. Ten or more classify a TaRca2 proteotype as RCA2\u03b1.",
   "identity_to_TaRca1_pct":"Highest percentage identity, over the aligned region, to any Chinese Spring proteotype of the TaRca1 loci of Caruana et al. (2022).",
   "identity_to_TaRca2_pct":"Highest percentage identity, over the aligned region, to any Chinese Spring proteotype of the TaRca2 loci of Caruana et al. (2022).",
   "identity_to_RCA2beta_pct":"Percentage identity to the canonical RCA2\u03b2 (RCA_PT006) over the aligned region; about 46 percent for the group 2 RCA-like lineage.",
   "assignment_evidence":"How the gene was assigned: directly, where a member protein model derives from a Chinese Spring anchor locus (TaRca1: TraesCS4A02G177600, TraesCS4B02G140200, TraesCS4D02G134900; TaRca2: TraesCS4A02G177500, TraesCS4B02G140300, TraesCS4D02G135000; Caruana et al. 2022); otherwise from sequence identity to those Chinese Spring proteotypes; or from homoeologous group 2 for the RCA-like lineage."}),
}

BY_BASE = {base: (title, desc) for base, title, desc in TABLES.values()}
FMT = {**numbering.SLABEL, **{f"T_{k}": v for k, v in numbering.TLABEL.items()}}
ORDER = [(i, key, title) for i, (key, _src, title) in enumerate(numbering.SUPP_TABLES, 1)]

def find(base, num):
    """Locate the source CSV for table <num>."""
    for pat in (f"TableS{num}_{base}.csv", f"TableS{num}_{base}*.csv", f"{base}.csv", f"{base}*.csv"):
        c = sorted(glob.glob(os.path.join(SRC, pat)))
        c = [p for p in c if os.path.isfile(p)]
        if c:
            return c[0]
    return None

SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
ok = True
for num, base, title in ORDER:
    if base not in BY_BASE:
        print(f"  NODESC S{num:<2} {base}: add its column descriptions to TABLES"); ok = False; continue
    desc = {c: d.format(**FMT) for c, d in BY_BASE[base][1].items()}
    src = find(base, num)
    if not src:
        print(f"  MISSING S{num:<2} {base}"); ok = False; continue
    with open(src, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    hdr, data = rows[0], rows[1:]
    missing = [c for c in hdr if c not in desc]
    extra   = [c for c in desc if c not in hdr]
    if missing or extra:
        print(f"  MISMATCH S{num}: undescribed={missing} not-in-file={extra}"); ok = False; continue
    dst = os.path.join(OUT, f"TableS{num}_{base}.csv")
    with open(dst, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([CITE])
        w.writerow([f"Supplemental Table S{num}. {title}"])
        w.writerow(hdr)
        w.writerows(data)
        w.writerow([])
        for c in hdr:
            w.writerow([c, desc[c]])
    print(f"  OK   S{num:<2} {base:<26} {len(data):>7} rows, {len(hdr):>2} cols -> {os.path.basename(dst)}")
sys.exit(0 if ok else 1)
