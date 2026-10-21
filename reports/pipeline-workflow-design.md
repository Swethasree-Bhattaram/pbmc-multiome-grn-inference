# Multimodal GRN pipeline: step contracts, gates, and method shortlist

Status: design document. No new runs are reported here; every number below is taken from
the existing reports in this repo, with the source file named.

Purpose: this repo has been iterated by changing inputs (reference choice, dataset size,
regulator list, target list) and reading the resulting GRN score. That loop cannot
converge, because the GRN score is dominated by choices made in the last two steps. This
document fixes the four steps, states what each step must produce before the next one may
run, states the gate that decides pass/fail, and lists recent methods to look at for each
step.

---

## 1. Diagnosis: what the existing runs say

### 1.1 The knob that has been varied does not matter

Same 10,844-cell all-cells matrix, TRRUST regulators, TRRUST targets; the only change is
which real RNA fed the imputation (`reports/comparison_report_trrust.md`):

| Experiment | Imputation reference | Recovered / 8,751 | Prec@100 | Prec@1000 |
|---|---|---|---|---|
| expA | SCEMENT-combined 3k + sub3k (5,422 cells) | 1603 | 0.0200 | 0.0100 |
| expB1 | 3k RNA only (2,711 cells) | 1657 | 0.0200 | 0.0100 |
| expB2 | sub3k RNA only (2,711 cells) | 1652 | 0.0100 | 0.0100 |

Reference choice moves recall by ±1.5% (1603 / 1657 / 1652). The full-10k variant repeats
the result (1957 / 1908 / 1939, `reports/full-10k/comparison_report_trrust_10k.md`). The
variable that has been under test is not a lever.

### 1.2 The multimodal path scores worse than the RNA-only baseline

| Run | Cells | Recall vs PBMC-TRRUST | Prec@500 | Prec@1000 |
|---|---|---|---|---|
| RNA-only baseline, no integration, no imputation (`reports/rna-baseline/REPORT_tfonly.md`) | 7,051 | 2655/8751 = 30.3% | 0.0380 | 0.0310 |
| expA, integrated + imputed (`reports/expA_REPORT.md`) | 10,844 | 1690/8751 = 19.3% | 0.0040 | 0.0060 |

Caveat, stated explicitly: this is **not** a controlled comparison. The baseline runs on
the 21,932-gene 3k+4k symbol intersection with the `tf_only` regulator list; expA runs on
the 36,601-gene axis with TRRUST regulators. The direction is consistent with the imputed
block diluting signal, but the controlled version of this comparison does not exist yet
(see §5, item 2).

### 1.3 The network formulation is broken independent of the biology

- 905,890 to 1,982,020 edges are emitted for 2,827 targets.
- Prec@1000 = 0.006–0.010: 6 to 10 of the top 1,000 edges are in the ground truth.
- Top edges (`reports/expA_REPORT.md`): RPL3→EEF1A1, B2M→MALAT1, HLA-B→MT-CYB,
  HLA-B→B2M, TPT1→EEF1A1, plus mutual pairs in both directions (CD8A↔CD8B, B2M↔RPL3,
  ANGPT2↔MCPH1). That is the ribosomal / mitochondrial / housekeeping axis plus a
  correlation-inflation artifact, not a regulatory program.

This is the documented behavior of imputation before network inference: imputation
inflates gene–gene correlations, which changes the reconstructed network structure and can
decrease reconstruction performance (Dijkstra et al., *Patterns* 2022,
[10.1016/j.patter.2021.100414](https://doi.org/10.1016/j.patter.2021.100414)). kNN
averaging is the most aggressive case, because a gene becomes nonzero in a cell if it is
nonzero in any one of its 20 neighbours.

### 1.4 Root cause: there is no measurement between steps

`grep -E "pearson|spearman|corrcoef|holdout|held_out|fidelity|silhouette|lisi"` over all
`*.py` in this repo returns **zero hits**. There is no imputation-fidelity check, no
batch-mixing check, no label-transfer check. The only signal produced anywhere is the
final GRN score — which §1.1 shows is insensitive to the upstream choices and §1.3 shows
is not a usable quality measure. Every upstream decision is therefore unfalsifiable, and
"try a different input" is the only available action. That is the problem this document
fixes.

---

## 2. The four steps and their contracts

Step 0 is not one of the four; it is a prerequisite, because without it none of the other
gates are reportable.

### Step 0 — Evaluation harness

**In:** any all-cells matrix plus its gene axis and cell mask.

**Out:**
- one script that accepts a matrix + gene axis + regulator list + target rule, and always
  also evaluates the RNA-only baseline through the identical harness;
- a fixed metric block: recall, prec@K for K ∈ {100, 500, 1000}, against both a whole-blood
  and a cell-type-matched ground truth.

**Gate:** the baseline reproduces 2655/8751 within seed noise on the same gene axis and
regulator list as the new runs. No downstream result is reportable until this exists.

**Why first:** it is the only way a later change can be called an improvement.

### Step 1 — Integration

**In:** per-dataset PCA (50 PCs), modality label, dataset ID. The paired 3k/10k multiome
cells are the only source of truth available and must be retained as a held-out set.

**Out (contract):**
- `H`: float32, `N_all × d`, `d ≥ 30` (scSAGA produces 30). Row order is recorded in a JSON
  map `dataset → [start, end)` and must not be re-derived later.
- per-dataset block slices, saved once.
- **no** gene axis, **no** reference choice, **no** imputation inside this artifact.

**Gate to enter step 2 (all three, on held-out paired cells):**

| # | Check | Definition | Pass criterion |
|---|---|---|---|
| a | Cross-modal kNN concordance | for a paired cell, take the 10 nearest neighbours in `H` of its ATAC half; fraction that are also within the 10 nearest of its RNA half | far above chance, and stable when the anchor is swapped |
| b | Label transfer | RNA-derived cell-type label assigned to the ATAC half of the same cell | > 0.85 on major PBMC lineages |
| c | Mixing vs separation | iLISI / cLISI (or kBET) per dataset | improves without (a) degrading |

Check (a) is anchor-free and is the direct measurement of whether `H` is usable as a
lookup space. It has never been computed here. Check (c) guards the opposite failure —
distribution matching that over-aligns biologically distinct states, which is the failure
mode the scHPGT paper (2026, §4.1) is written against.

Do **not** use scSAGA's "global alignment score" as a gate: on 5 distinct 10x experiments
it moved 0.0842 → 0.1700 (2.02×) from the anchor alone, with no bug involved.

**On failure:** change the integrator. Do not touch the GRN step.

### Step 2 — Imputation

**In:** `H`, the step-1 block map, and the RNA reference on one declared normalization.

**Out (contract):**
- `X`: float32, `N_all × G`, with `G` = the intersection gene axis, written once to
  `genes.txt` and never re-derived.
- a boolean column mask per cell: real vs imputed.
- one declared normalization, **identical for real and imputed rows**. The current imputed
  block (mean 0.0911, nonzero fraction 0.490, `reports/expA_REPORT.md`) has not been shown
  to match the real log1p(CPM/1e4) rows, and any mismatch becomes a batch covariate for
  GRNBoost2.
- the imputer's fitted reference set, frozen, so the held-out test is honest.

**Gate to enter step 3:**

| # | Check | Definition | Pass criterion |
|---|---|---|---|
| a | ATAC→RNA fidelity | hold out paired 3k (and 10k) cells, impute from ATAC only, per-gene Pearson + Spearman on held-out cells, plus pseudobulk correlation | at or above the published ceiling for this task; ≫ 0 |
| b | Densification | nonzero fraction of imputed rows vs real rows | not materially above real |
| c | Cross-modal consistency | imputed expression vs gene activity of linked peaks, same cell | positive and gene-specific |

Calibrate (a) against Polarbear's published numbers: ATAC→RNA gene-wise Pearson ≈ 0.22
(Polarbear) vs 0.155 (BABEL) on SNARE-seq. The ceiling for this task is low in absolute
terms; a run reporting ~0.05 is producing noise, and step 3 on that matrix is uninterpretable.

**On failure:** replace the imputer (§4.2) or drop step 2 and give step 3 gene activity
alongside real RNA. Do not run step 3 on a matrix that failed its own gate.

### Step 3 — GRN inference

**In:** `X` plus the column mask, and — strongly recommended — gene activity from ATAC as
a second matrix.

**Out (contract):**
- edge list `TF, target, importance`, plus `cells_used` and the exact software stack.
- **three seeds**, not one, so stability is measurable.
- networks computed **per cell type or on pseudobulk**, not one network over 10,844
  heterogeneous cells.
- regulators = a curated TF list; targets = a stated rule with a stated count.

**Gate to enter step 4:**

| # | Check | Pass criterion |
|---|---|---|
| a | Seed stability | Jaccard of top-10k edges across 3 seeds reported; differences below that noise floor are not interpreted |
| b | Lineage sanity | PAX5/EBF1→MS4A1/CD79A in B; SPI1/CEBPA→LYZ/S100A8 in myeloid; TBX21/GNK in NK/T; CD3D/CD3E in T — present in the *cell-type* network's top edges |
| c | Beats baseline | beats the RNA-only baseline on the step-0 harness by a stated margin |

Seed sensitivity is not hypothetical: GRETA found multimodal GRNs "highly sensitive to
method choices, such as changes in random seeds, or replacing steps in the inference
pipelines" (§4.4).

**On failure:** go imputation-free (LINGER / SCENIC+ / HuMMuS, §4.3) rather than tuning
GRNBoost2.

### Step 4 — Evaluation and validation

**In:** per-cell-type edge lists from step 3.

**Out (contract):** a metric block that never rests on TRRUST edge overlap alone.

| Axis | Measurement |
|---|---|
| Causal | predicted response to TF perturbation vs held-out Perturb-seq (scPerturb / Norman / Replogle) |
| Cell-type-resolved | overlap with cell-type-specific ChIP / ReMap / ChIP-Atlas edges, not the 96,846-edge whole-blood blob |
| Regulon-level | SCENIC-style AUCell activity of top regulons per cell type |
| Predictive | does the network rank a held-out TF's true targets above background |

**Gate:** prec@500 > 0.10 against cell-type-matched ground truth. Today's 0.004–0.038 means
the top of the ranking carries no information and no biological claim can be made from it.

The 96,846-edge "PBMC-Blood" truth set is part of why precision reads as zero; it is a
whole-blood resource being scored against a sorted-PBMC network.

---

## 3. Consequences for how work is sequenced

1. Build step 0 first.
2. Develop steps 1 and 2 entirely against their gates — no GRNBoost2, no dask, days not
   weeks of compute.
3. Run the GRN once the matrix has passed its gate, with 3 seeds, per cell type.
4. Only then interpret biology, and only if step 4's gate passes.

This inverts the current loop, where the GRN was re-run to test upstream changes.

---

## 4. Papers, by step

Metadata checked against Crossref / Europe PMC. Each entry names the step it maps to and
the specific reason to read it.

### 4.1 Step 1 — integration

1. **Liu et al., "Multitask benchmarking of single-cell multimodal omics integration
   methods", *Nature Methods* 2025.**
   [10.1038/s41592-025-02856-3](https://doi.org/10.1038/s41592-025-02856-3)
   Read first. Systematic categorization plus benchmarking of current methods across
   dimension reduction, batch correction, cell-type classification, clustering,
   **imputation**, feature selection and spatial registration, with a decision guideline.
   It gives the selection rule for both step 1 and step 2 instead of guessing. Companion
   benchmark: "Benchmarking single-cell multi-modal data integrations", *Nature Methods*
   2025, [10.1038/s41592-025-02737-9](https://doi.org/10.1038/s41592-025-02737-9)
   (40 methods over paired, unpaired and mosaic data).

2. **Qi et al., "scGT: integration algorithm for single-cell RNA-seq and ATAC-seq based on
   graph transformer", *Bioinformatics* 2025.**
   [10.1093/bioinformatics/btaf357](https://doi.org/10.1093/bioinformatics/btaf357)
   Graph transformer over correlation-strengthened graphs; paired and unpaired; evaluated
   on label transfer and on millions of cells. Maps directly onto step-1 gates (a) and (b).

3. **"SuperMap: bridging unpaired single-cell multimodal data for integrative analyses",
   *PNAS* 2026.**
   [10.1073/pnas.2505182123](https://doi.org/10.1073/pnas.2505182123)
   Learns cross-modal mappings directly from unpaired data. Becomes relevant the moment
   ATAC-only datasets are added to the integration.

4. **"scHPGT: regulatory-prior-guided attention preserves biological structure during
   unpaired single-cell RNA–ATAC integration", *Bioinformatics* 2026.**
   [10.1093/bioinformatics/btag696](https://doi.org/10.1093/bioinformatics/btag696)
   Argues explicitly against framing integration as distribution matching, which can
   over-align condition-specific and modality-specific states. Use its critique as the
   justification for gate (c).

If a single model should own steps 1 and 2 together: **MultiVI** (Ashuach, Gabitto et al.,
*Nature Methods* 2023, [10.1038/s41592-023-01909-9](https://doi.org/10.1038/s41592-023-01909-9))
jointly models RNA + ATAC and infers the missing modality probabilistically.

### 4.2 Step 2 — imputation

1. **"scGALA: advancing graph link prediction-based cell alignment for comprehensive data
   integration and harmonization", *Nature Communications* 2025.**
   [10.1038/s41467-025-66644-5](https://doi.org/10.1038/s41467-025-66644-5)
   Closest drop-in for the current design: an alignment step plus a graph-attention
   `GATMapper` that generates RNA for ATAC-only cells from the learned alignment. Same
   two-stage structure already in this repo, with a defensible mapper in place of
   `exp(-distance)` neighbour averaging.

2. **"scDCT: a conditional diffusion-based deep learning model for high-fidelity
   single-cell cross-modality translation", *Briefings in Bioinformatics* 2025.**
   [10.1093/bib/bbaf400](https://doi.org/10.1093/bib/bbaf400)
   Bidirectional RNA↔ATAC with explicit uncertainty; benchmarked against BABEL, Polarbear,
   JAMIE and scButterfly, with better preservation of cell-subpopulation topology.

3. **"Cisformer: a scalable cross-modality generation framework for decoding
   transcriptional regulation at single-cell resolution", *Genome Biology* 2025.**
   [10.1186/s13059-025-03823-z](https://doi.org/10.1186/s13059-025-03823-z)
   10,000 active CREs per cell, released pretrained weights. The option to take if the
   pipeline scales past ~30k cells.

4. **"scButterfly: a versatile single-cell cross-modality translation method via
   dual-aligned variational autoencoders", *Nature Communications* 2024.**
   [10.1038/s41467-024-47418-x](https://doi.org/10.1038/s41467-024-47418-x)
   Simplest method that handles paired and unpaired settings. Useful as a control arm
   precisely because it is not a diffusion model.

Read before building either of the above:

- **Polarbear** (Zhang, Meng-Papaxanthos, Vert, Noble, "Multimodal single-cell translation
  and alignment with semi-supervised learning", *J Comput Biol* 2022,
  [10.1089/cmb.2022.0264](https://doi.org/10.1089/cmb.2022.0264)): publishes ATAC→RNA
  gene-wise Pearson ≈ 0.22 (Polarbear) vs 0.155 (BABEL) on SNARE-seq — the accuracy
  ceiling, in absolute terms, is low.
- **Dijkstra et al., "Effect of imputation on gene network reconstruction from
  single-cell RNA-seq data", *Patterns* 2022**
  ([10.1016/j.patter.2021.100414](https://doi.org/10.1016/j.patter.2021.100414)):
  imputation inflates gene–gene correlations, alters predicted network structure and can
  decrease reconstruction performance. This is the mechanism behind §1.3.

### 4.3 Step 3 — GRN inference, and the imputation-free alternative

1. **Yuan et al., "Inferring gene regulatory networks from single-cell multiome data using
   atlas-scale external data" (LINGER), *Nature Biotechnology* 2025.**
   [10.1038/s41587-024-02182-7](https://doi.org/10.1038/s41587-024-02182-7)
   Paired RNA + ATAC, TF-motif prior, atlas-scale bulk data as manifold regularization;
   reports a 4- to 7-fold relative accuracy increase over existing methods. This is the
   paper that makes step 2 optional: keep chromatin in the model instead of translating it
   into expression first.

2. **Bravo Gonzalez-Blas et al., "SCENIC+: single-cell multiomic inference of enhancers and
   gene regulatory networks", *Nature Methods* 2023.**
   [10.1038/s41592-023-01938-4](https://doi.org/10.1038/s41592-023-01938-4)
   Still the reference multiome GRN method, and it returns **regulons** — a far more stable
   evaluation unit than a 2M-edge ranked list.

3. **"HuMMuS: molecular mechanisms reconstruction from single-cell multi-omics data",
   *Bioinformatics* 2024.**
   [10.1093/bioinformatics/btae143](https://doi.org/10.1093/bioinformatics/btae143)
   Multilayer network fusing RNA, ATAC and PPI, with no imputation step in the path.

4. **"scMEGA: single-cell multi-omic enhancer-based gene regulatory network inference",
   *Bioinformatics Advances* 2023.**
   [10.1093/bioadv/vbad003](https://doi.org/10.1093/bioadv/vbad003)
   Lightweight Bioconductor route (integration → peak–gene links → GRN). The fast,
   inspectable option when a deep model is not warranted.

Context: geneRNIB (§4.4, ref 2) reports that simple models with fewer assumptions
outperform complex ones on this task. GRNBoost2 is not the first thing to replace — the
input matrix and the evaluation are.

### 4.4 Step 4 — evaluation and validation

1. **Badia-i-Mompel et al., "Comparison and evaluation of methods to infer gene regulatory
   networks from multimodal single-cell data" (GRETA), bioRxiv 2025.**
   [10.1101/2024.12.20.629764](https://doi.org/10.1101/2024.12.20.629764)
   The single most relevant paper here. GRNs from multimodal data are highly sensitive to
   method choices including random seeds and step substitution; they perform well on
   predictive tasks, only partially recover known interactions, and "struggled to capture
   causal relationships from perturbation assays". Shipped as a modular Snakemake pipeline,
   so each of the four steps becomes a swappable module behind a fixed harness.

2. **Nourisa et al., "geneRNIB: a living benchmark for gene regulatory network inference",
   bioRxiv 2025.**
   [10.1101/2025.02.25.640181](https://doi.org/10.1101/2025.02.25.640181)
   12 methods over 11 datasets and thousands of perturbation scenarios; introduces
   context-specific metrics and evaluation in the absence of a true reference network. The
   methodology for step 4's "causal" and "predictive" axes.

3. **"A large-scale benchmark for network inference from single-cell perturbation data",
   *Communications Biology* 2025.**
   [10.1038/s42003-025-07764-y](https://doi.org/10.1038/s42003-025-07764-y)
   Real perturbation-based evaluation, showing that synthetic benchmarks do not reflect
   real-world performance.

4. **Peidli et al., "scPerturb: harmonized single-cell perturbation data", *Nature Methods*
   2024.**
   [10.1038/s41592-023-02144-y](https://doi.org/10.1038/s41592-023-02144-y)
   The resource for causal ground truth. Use this instead of adding further TRRUST
   edge-overlap runs.

Also worth reading, on the evaluation unit itself:
"Population-level comparisons of gene regulatory networks modeled on high-throughput
single-cell transcriptomics data", *Nature Computational Science* 2024
([10.1038/s43588-024-00597-5](https://doi.org/10.1038/s43588-024-00597-5)) — argues that
single-cell-level GRN inference is not reproducible and motivates population/pseudobulk
comparison. This is the paper behind "compute the network per cell type on pseudobulk".

---

## 5. Minimum changes, in order

1. **Add the step-2 held-out fidelity check** on the paired 3k cells. Nothing in the repo
   currently measures imputation quality; the paired multiome makes this test cheap and
   ground-truthed. This is the single highest-value addition.
2. **Move the RNA-only baseline into the same harness** — same gene axis, target list and
   regulator list as the multimodal runs — and report the delta. Today's 30.3% vs 19.3%
   is suggestive but uncontrolled (§1.2).
3. **Stop reporting prec@K on a 1–2M edge network.** Per-cell-type networks,
   top-N-per-target, cell-type-matched ground truth.
4. **Run 3 seeds before believing any network difference.** GRETA shows the seed alone
   moves the network more than most method changes.

---

## 6. Cited sources

| Step | Reference | Identifier |
|---|---|---|
| 1 | Liu et al., *Nat Methods* 2025, multitask benchmarking | 10.1038/s41592-025-02856-3 |
| 1 | Benchmarking single-cell multi-modal data integrations, *Nat Methods* 2025 | 10.1038/s41592-025-02737-9 |
| 1 | Qi et al., scGT, *Bioinformatics* 2025 | 10.1093/bioinformatics/btaf357 |
| 1 | SuperMap, *PNAS* 2026 | 10.1073/pnas.2505182123 |
| 1 | scHPGT, *Bioinformatics* 2026 | 10.1093/bioinformatics/btag696 |
| 1–2 | Ashuach, Gabitto et al., MultiVI, *Nat Methods* 2023 | 10.1038/s41592-023-01909-9 |
| 2 | scGALA, *Nat Commun* 2025 | 10.1038/s41467-025-66644-5 |
| 2 | scDCT, *Brief Bioinform* 2025 | 10.1093/bib/bbaf400 |
| 2 | Cisformer, *Genome Biol* 2025 | 10.1186/s13059-025-03823-z |
| 2 | scButterfly, *Nat Commun* 2024 | 10.1038/s41467-024-47418-x |
| 2 | Zhang et al., Polarbear, *J Comput Biol* 2022 | 10.1089/cmb.2022.0264 |
| 2–3 | Dijkstra et al., *Patterns* 2022, imputation and network reconstruction | 10.1016/j.patter.2021.100414 |
| 3 | Yuan et al., LINGER, *Nat Biotechnol* 2025 | 10.1038/s41587-024-02182-7 |
| 3 | Bravo Gonzalez-Blas et al., SCENIC+, *Nat Methods* 2023 | 10.1038/s41592-023-01938-4 |
| 3 | HuMMuS, *Bioinformatics* 2024 | 10.1093/bioinformatics/btae143 |
| 3 | scMEGA, *Bioinformatics Advances* 2023 | 10.1093/bioadv/vbad003 |
| 4 | Badia-i-Mompel et al., GRETA, bioRxiv 2025 | 10.1101/2024.12.20.629764 |
| 4 | Nourisa et al., geneRNIB, bioRxiv 2025 | 10.1101/2025.02.25.640181 |
| 4 | Large-scale benchmark, single-cell perturbation, *Commun Biol* 2025 | 10.1038/s42003-025-07764-y |
| 4 | Peidli et al., scPerturb, *Nat Methods* 2024 | 10.1038/s41592-023-02144-y |
| 4 | Population-level comparisons of GRNs, *Nat Comput Sci* 2024 | 10.1038/s43588-024-00597-5 |
