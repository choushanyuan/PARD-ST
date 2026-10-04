# PARD-ST

Analysis code for the spatial transcriptomics (10x Visium, FFPE) study of prostate
cancer — **PARD-ST**.

Seven patient samples were profiled: `P1`–`P5` (sample IDs `PR_P1`–`PR_P5`),
`P6` (sample ID `P60`) and `P7` (sample ID `P616`). Each sample is processed on its
own and the sections are then integrated prior to the cross-patient analysis.

---

## Repository layout

| Path | Content |
| --- | --- |
| `P1/` … `P7/` | per-patient analysis (one directory per sample) |
| `integration/` | SCTransform integration of the seven sections, GSVA heatmaps, spatial-continuity / transcriptome-diversity metrics |
| `survival_analysis/` | bulk-cohort (Changhai and TCGA) survival analysis of the pseudotime-derived signatures |
| `survival_analysis/tools/` | small helper scripts (signature table export, PDF splitting) |
| `survival_analysis/code_backup/` | earlier revisions of `01`–`04`, kept for reference |

### Per-patient files (identical layout in every `P1/` … `P7/`)

| File | Step |
| --- | --- |
| `spaceranger_notes.md` | Space Ranger 2.0 alignment (`spaceranger count`, probe set, Loupe alignment) |
| `seurat_clustering_annotation.md` | object creation, SPROD de-noising, QC, clustering, pathologist annotation, scRNA-seq label transfer |
| `estimate_notes.md` | ESTIMATE tumour purity / stromal / immune scores (`P2` and `P7` have no estimate step) |
| `infercnv_notes.md` | inferCNV input preparation, run, dendrogram cutting and CNV subclones |
| `step1_pseudotime.R` | monocle pseudotime (`DDRTree`) and projection of pseudotime/state onto the sections |
| `step2_transcription_factors.R` | DoRothEA regulon + VIPER transcription-factor activity, clustered heatmap |
| `step3_pie_chart.R` | composition of the monocle states across the annotated clusters (pie charts) |
| `step4_differential_genes.R` | `FindMarkers` between states and volcano plots |
| `P?_cellchat.R` | CellChat cell–cell communication (spatial mode) and pathway visualisation |

The scripts of a sample are meant to be run in the order above: `step2`–`step4`
reuse the `ST_obj` object created by `step1` (they do not load it themselves), and
`step3` reads `Cellratio.RData` written by `step1`.

---

## Integration (`integration/`)

| File | Content |
| --- | --- |
| `prad_intergrated.R` | SCTransform anchor integration of all seven sections, PCA / clustering / UMAP, `FindAllMarkers` on `integrated_snn_res.0.2` |
| `heatmap_by_patient.R` | GSVA of hallmark and cell-junction gene sets, averaged per tumour subtype and per patient |
| `cell_junction_pathway_heatmap.R` | GSVA of cell-junction pathways across the tumour subtypes |
| `spatial_continuity_and_transcriptome_diversity.md` | per-spot transcriptome diversity (MAD of the Pearson correlation) and spatial continuity, plus the scatter plot |

Tumour subtypes used throughout the integrated analysis:

```
subtype_low    = integrated_snn_res.0.4 cluster 1
subtype_middle = integrated_snn_res.0.4 clusters 0, 2
subtype_high   = integrated_snn_res.0.4 clusters 9, 12
```

## Survival analysis (`survival_analysis/`)

Run in numeric order; each script reads what the previous one wrote to
`<project_root>/result/` and writes figures to `<project_root>/figure/`.
`04b` and `07` derive the project root as the **parent of their own directory**, so
`result/` and `figure/` are expected one level above the scripts.

| File | Content |
| --- | --- |
| `01_calculate_scores.R` | ssGSEA scores (`Overlap_genes`, `Up_score`, `Down_score`) and ESTIMATE scores for the Changhai and TCGA cohorts, for both FPKM and TPM inputs |
| `02_kaplan_meier.R` | Kaplan-Meier curves for HSPD1, FASN, PKP1, `Up_score`, `Down_score` under five cutoffs (Median, Q3_Upper, Q1_Lower, Q1_vs_Q3, Optimal) |
| `03_multivariable_cox.R` | multivariable Cox model with age stratification at 65 years, forest plots and model diagnostics (proportional hazards, linearity) |
| `04_gene_contributions.R` | per-gene univariate Cox (continuous expression) forest plots |
| `04b_gene_contributions_median.R` | the same per-gene analysis with a median high/low split |
| `05_investigate_variables.R` | diagnostic dump of the clinical columns (TNM stage, surgical margin) and the recoding strategy |
| `06_compare_cindex.R` | C-index of `Clinical only` vs `Clinical + Score`, likelihood-ratio test and comparison bar charts |
| `07_raincloud_CRISP3_AMACR_MYC.R` | raincloud plots of CRISP3 / AMACR / MYC expression across three Gleason strata |

Signature direction (as documented in `tools/export_signature_gene_tables.R`):

```
Up signature   = down_genes  (n = 95)   # pseudotime dedifferentiation up-regulated genes
Down signature = up_genes    (n = 246)
```

The names inside the `.RData` files follow the raw pseudotime direction and are
therefore inverted with respect to the manuscript definition — see the comments at
the top of `tools/export_signature_gene_tables.R` before changing them.

---

## Requirements

* **R** (≥ 4.1) with `Seurat`, `monocle`, `infercnv`, `estimate`, `dorothea`, `viper`,
  `CellChat`, `GSVA`, `msigdbr`, `fgsea`, `survival`, `survminer`, `forestploter`,
  `ggdist`, `ggpubr`, `patchwork`, `pheatmap`, `tidyverse`, `showtext`
* **Space Ranger 2.0** and **Loupe Browser** for the upstream alignment and spot selection
* **Python 3** with `PyMuPDF` (`fitz`) and `pikepdf` — only for
  `survival_analysis/tools/split_pdf_for_illustrator.py`

