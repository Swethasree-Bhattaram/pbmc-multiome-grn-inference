# Multimodal GRN pipeline: a workflow with a measurement at every step

This document describes four steps — integration, imputation, GRN inference, validation —
written so that each step can be built, tested and rejected on its own evidence.

Two rules run through the whole document:

1. Every step has one input, one output, and one measurement.
2. A step must pass its measurement before the next step is run.

Rule 2 is not tidiness. In a four-step pipeline the only number produced at the end is the
network score, and that score cannot say which of the four steps failed. A mistake in step 1
does not disappear; it gets organised by step 2 and ranked confidently by step 3. Measuring
each step is the only way to know where a problem is.

Every step section below has the same five parts, in the same order:

| Part | Question it answers |
|---|---|
| What it does | one paragraph, no jargon |
| Input | what must exist before the step can start |
| Output | what exists afterwards, and what it must never contain |
| Methods | what can do it, with citations |
| How we measure it | the metrics, the ground truth they need, and the tool |
| Gate | the pass condition, and what to do on failure |

---

## 0. Vocabulary

Ten terms are used throughout. Nothing else is assumed.

| Term | Plain meaning |
|---|---|
| RNA (scRNA-seq) | which genes each cell is transcribing |
| ATAC (scATAC-seq) | which DNA is open — what *could* be transcribed |
| multiome | RNA and ATAC measured in the same cell |
| modality | one kind of measurement (RNA is one, ATAC is another) |
| embedding | ~30 numbers per cell, so cells can be compared by distance |
| integration | one shared embedding for cells measured different ways |
| imputation | filling in a measurement a cell does not have |
| gene activity | a per-gene score computed from open chromatin near that gene |
| GRN | a table of TF → target gene links, with a score |
| regulon | a TF plus the genes it drives, scored per cell |

Method names later in the document are only different ways of doing one of these ten things.

---

## 1. The workflow at a glance

```
1. INTEGRATION             2. IMPUTATION              3. GRN                 4. VALIDATION
   IN                         IN                         IN                     IN
   · one PCA per dataset      · H  (from step 1)         · the matrix X         · edge tables
   · which modality a cell is · measured RNA: the         · TF list              · independent
     (RNA / ATAC)               values being copied      · target rule            evidence
   · which dataset it is        onto ATAC cells          · cells grouped        · a negative
   · paired cells set aside   · which cells need           by cell type           control
     as the answer key          values (ATAC-only)
   OUT                        OUT                        OUT                    OUT
   · H: cells × 30 numbers    · X: cells × all genes     · TF, target, score    · AUPRC, early
   · dataset → row range      · mask: measured/guessed   · cells used, 3 seeds    precision,
   · no genes, no imputation  · one normalization        · one network            stability
                                                          per cell type         · a verdict
```

Five phrases in step 2 are worth defining once — three in its input, two in its output:

| Phrase in the diagram | What it actually is |
|---|---|
| **measured RNA** (also "real RNA") | the gene-expression matrix of the cells that *have* RNA — the multiome cells. It is the only source of the values that get copied onto ATAC cells. Called "measured" only to distinguish it from the guessed values. |
| **which cells need values** ("ATAC row list") | the row numbers of `X` that belong to ATAC-only cells. They have no measurement, so they are the rows to be filled in. Without this list you cannot say which rows of the output are guesses. |
| **paired cells set aside as the answer key** ("held out", "frozen") | a set of cells measured *both* ways, whose RNA we hide from the imputer and compare against afterwards. "Frozen" means: chosen once, recorded, and never allowed into the imputer's reference set or any fitting — otherwise the comparison is against something the method already saw. This is the test set for steps 1 and 2. |
| **mask: measured / guessed** | a true/false grid the same shape as `X`. True where the number came from a real measurement, false where it is a guess. It is how every later step — and every metric — can ask "was this value measured, or invented?", and it is what lets step 3 run with and without the guessed rows. |
| **one normalization** | exactly one declared transformation from raw counts to comparable numbers (e.g. `log1p(CPM/1e4)`), written to `norm.txt`, applied so that measured rows and guessed rows sit on the **same scale**. If the imputer emits something else — a probability in [0, 1], say — it must be mapped into that space. Otherwise the two blocks are separable, and the network model learns "which block is this" instead of biology. |

Each column's output is literally the next column's input. That is the whole design: there
are only two objects to inspect, `H` and `X`, and both can be inspected without running a
network.

### What physically crosses each boundary

| After | File | Shape | Must never contain |
|---|---|---|---|
| 1 | `H.npy` | cells × 30, float32 | genes, peaks, any imputation |
| 1 | `blocks.json` | dataset → row range | — |
| 2 | `X.npy` | cells × genes, float32 | a second gene axis, mixed scales |
| 2 | `mask.npy` | cells × genes, bool | — |
| 2 | `genes.txt`, `norm.txt` | one gene axis, one normalization | — |
| 3 | `edges.csv` | TF, target, importance | unrecorded software or seed |
| 3 | `run_info.txt` | cells, columns, seeds, versions | — |
| 4 | `metrics.json` | same metrics for every method | per-method metric choices |

---

## 2. Ground truth: the four kinds, and what each is good for

Every metric in this document needs one of these. Two are free with multiome data, which is
why steps 1 and 2 can be graded strictly.

| # | Kind | What it is | Good for | Limitation |
|---|---|---|---|---|
| 1 | Paired cells | the same cell measured both ways | steps 1 and 2 | only available for multiome |
| 2 | Held-out labels | known cell types, known batches | step 1 | labels can be wrong at boundaries |
| 3 | External prior | ChIP-seq, curated databases, Perturb-seq | step 3 | incomplete, not cell-type specific |
| 4 | Simulation | a network we specify ourselves | the whole pipeline | artificial |

Simulation is not a substitute for biology. It is the only place where we can ask whether
the pipeline recovers a network that is definitely there. That question is answered before
the pipeline is pointed at real data (section 6).

---

## 3. Step 1 — Integration

### What it does

Two measurements of the same cells must end up on one shared map. An ATAC cell and the RNA
cell from the same B cell should be close together; two different cell types should not be.

### Input

| Item | Shape / form | Why it is needed |
|---|---|---|
| per-dataset PCA | cells × 50 | integration never sees genes or peaks, only 50 summaries per cell |
| modality label | RNA / ATAC | so the method knows which side a cell came from |
| dataset label | string per dataset | needed for the batch-mixing metrics |
| paired cells | a cell index, held out | the only ground truth for this step; must not be used to fit |

Note on ATAC feature space: the scIB benchmark found that ATAC integration performance is
strongly affected by the choice of feature space, more so than for RNA. Fix the ATAC
feature space (peaks, or gene activity) before comparing integrators, and record it.

### Output

| Item | Shape / form | Notes |
|---|---|---|
| `H.npy` | cells × 30, float32 | the only thing step 2 consumes |
| `blocks.json` | dataset → `[start, end)` | row order must never be re-derived |
| per-dataset slices | one file each | saved once |

Deliberately absent: no gene axis, no imputation, no decision about which RNA dataset feeds
the imputation. Those belong to step 2. Putting them here is how a pipeline loses the
ability to say which step caused a change.

### Methods

| Method | What it does | Choose it when | Reference |
|---|---|---|---|
| scSAGA | optimal-transport alignment to an anchor dataset | current baseline | Aluru lab, scSAGA |
| scGT | graph transformer over correlation-strengthened graphs; paired and unpaired; label transfer | scaling to millions of cells | Qi et al., *Bioinformatics* 2025, 10.1093/bioinformatics/btaf357 |
| scMODAL | aligns using a few known linked features, GAN-based | few feature links available | *Nat Commun* 2025, 10.1038/s41467-025-60333-z |
| SuperMap | learns cross-modal mappings directly from unpaired data | ATAC-only datasets are added | *PNAS* 2026, 10.1073/pnas.2505182123 |
| scHPGT | regulatory-prior-guided attention; resists over-alignment | cell states must stay distinct | *Bioinformatics* 2026, 10.1093/bioinformatics/btag696 |
| MultiVI | one probabilistic model for RNA + ATAC, infers the missing modality | steps 1 and 2 should be one model | Ashuach, Gabitto et al., *Nat Methods* 2023, 10.1038/s41592-023-01909-9 |

Read before choosing: Liu et al., *Nat Methods* 2025
(10.1038/s41592-025-02856-3) benchmarks 40 methods across these exact tasks — dimension
reduction, batch correction, clustering, label transfer and imputation — and publishes a
decision guideline. Companion: *Nat Methods* 2025, 10.1038/s41592-025-02737-9 (40 methods,
paired/unpaired/mosaic).

### How we measure it

Four metrics. The first needs paired cells and is the strictest; the others need labels.

| Metric | What it asks | Ground truth | Good looks like | Where |
|---|---|---|---|---|
| FOSCTTM | for each cell, what fraction of other cells are closer to it than its true cross-modal partner is? | paired cells (held out) | 0 = perfect, 0.5 ≈ random. Published on SNARE-seq: 0.14–0.16 (MMD-MA), 0.22 (UnionCom), 0.39 (no normalization) | ~10 lines of numpy |
| label conservation | are cell types still separated? ARI, NMI, cLISI, label silhouette, isolated-label F1, graph connectivity | cell-type labels | closer to unintegrated labels than to a shuffled control | `scib` |
| batch correction | are datasets mixed? kBET acceptance, iLISI, batch silhouette | dataset IDs | improves, while label metrics do not degrade | `scib` |
| label transfer | can a label learned on RNA be assigned to the ATAC half of the same cell? | paired cells + labels | > 0.85 on major lineages | own script |

The `scib` package implements 14 evaluation metrics covering the middle two rows. The scIB
benchmark itself (Luecken et al., *Nat Methods* 2022, 10.1038/s41592-021-01336-8) evaluated
68 method/preprocessing combinations on 85 batches and reported that some preprocessing
choices push methods to remove batch at the cost of destroying biological variation. Report
both halves separately for that reason; a single combined score can hide the trade.

FOSCTTM definition, for implementation: for each cell in one modality, rank all cells of the
other modality by distance in `H`, then compute the fraction of them that are closer than
the cell's true partner. Average over all cells and over both directions. Zero means every
cell's nearest neighbour across modalities is its true partner.

### Gate

Pass if all three hold:

1. FOSCTTM clearly below the shuffled-partner control (the same computation with the pairing
   permuted).
2. Label transfer above 0.85 on major lineages.
3. scIB batch metrics improve while scIB label metrics do not degrade.

Do not use scSAGA's "global alignment score" as a gate. It is computed relative to the
chosen anchor, so it measures the anchor as much as the embedding: on the same data it can
move by a factor of two when only the anchor changes, with no bug involved.

On failure: change the integrator. Do not touch later steps.

---

## 4. Step 2 — Imputation

### What it does

ATAC cells have no RNA measured. Step 2 produces a gene vector for them, so that every cell
— paired or not — can enter one matrix.

The simplest version: take one ATAC cell, find the 20 measured-RNA cells closest to it in
`H`, average their expression, and call that the cell's expression. That is the whole idea,
and it is also where the risk lives.

### Why averaging is dangerous

Averaging is smoothing, and smoothing manufactures correlation. Four consequences, in order:

| # | What happens | Why it matters |
|---|---|---|
| 1 | Zeros disappear | an average of 20 cells is rarely zero, so guessed cells look "richer" than real ones |
| 2 | Correlations inflate | any two genes present in the same neighbourhood become correlated — because we averaged them, not because a TF links them |
| 3 | The block is detectable | guessed rows have a different distribution from real rows; a model can separate them |
| 4 | The GRN chases 1–3 | a model that ranks genes by how well they predict each other ranks inflated correlations first |

The visible signature of consequence 4 is a network whose strongest links are housekeeping
genes — ribosomal, mitochondrial, HLA — because those are high and correlated in every
neighbourhood.

This is not a hypothetical concern. Dijkstra et al., *Patterns* 2022
(10.1016/j.patter.2021.100414) evaluated networks before and after imputation and found
inflated gene–gene correlations, altered network structure, and a possible decrease in
reconstruction performance.

Consequences 1–3 can all be measured directly. That is what the metrics below do.

### Input

| Item | Shape / form | Why it is needed |
|---|---|---|
| `H.npy` + `blocks.json` | cells × 30 | the lookup space, from step 1 |
| measured RNA | cells × genes, **one** gene axis | the values being copied onto the ATAC cells |
| one normalization | declared, e.g. log1p CPM/1e4 | measured and guessed rows must be on the same scale |
| which cells need values | row indices of the ATAC-only cells | these are the rows of `X` being filled in |
| answer-key cells | paired cells, fixed and excluded from fitting | the hidden RNA we compare against afterwards; if they enter the reference, the test is void |

Two details that cause silent failures. First, the gene axis is written to a file once and
never re-derived: two experiments on 36,601 and 21,932 genes produce numbers that cannot be
compared. Second, measured and guessed rows must share one normalization, because a mismatch
between the two becomes a batch covariate that the network model will happily learn.

A one-line check catches both: for each gene, compare mean, variance and nonzero fraction
between the measured rows and the guessed rows. They will not be identical — the guesses are
smoother — but a systematic shift across most genes means the two blocks are not on the same
scale, and step 3 will partly be learning the block, not the biology.

### Output

| Item | Shape / form | Notes |
|---|---|---|
| `X.npy` | cells × genes, float32 | the only input to step 3 |
| `genes.txt` | one line per gene | written once |
| `mask.npy` | cells × genes, bool | per cell per gene: measured or guessed |
| `norm.txt` | the one normalization | applied to both kinds of row |
| reference set actually used | cell indices | proof the answer key was never inside it |

The mask is what makes step 2 testable. Without it, nothing downstream can point at a row
and say "that one is an average", so nobody can ask whether the averages helped. With it,
step 3 can be run with and without the guessed rows — the comparison that decides whether
step 2 earns its place.

### Methods

| Method | In plain words | Choose it when | Reference |
|---|---|---|---|
| kNN average | average the 20 nearest RNA cells | simplest; see the table above for the cost | this repo |
| scGALA | step-1 alignment plus a learned mapper (graph attention) | closest drop-in for an alignment-then-map design | *Nat Commun* 2025, 10.1038/s41467-025-66644-5 |
| scDCT | conditional diffusion, translation both directions | topology preservation matters most | *Brief Bioinform* 2025, 10.1093/bib/bbaf400 |
| Cisformer | pretrained; 10,000 open regions per cell | more than ~30k cells | *Genome Biol* 2025, 10.1186/s13059-025-03823-z |
| scButterfly | two aligned VAEs, paired or unpaired | simplest learned model, good control arm | *Nat Commun* 2024, 10.1038/s41467-024-47418-x |
| MultiVI | one model does steps 1 and 2 | paired data only | *Nat Methods* 2023, 10.1038/s41592-023-01909-9 |
| Polarbear | the original semi-supervised translator | the number to beat | *J Comput Biol* 2022, 10.1089/cmb.2022.0264 |

The published ceiling for this task is low. Polarbear reports ATAC→RNA gene-wise Pearson of
about 0.22, BABEL about 0.155, on SNARE-seq. Two consequences: a method claiming 0.9
correlation is more likely leaking than excellent, and a method reporting 0.05 has produced
nothing. Both are knowable in an afternoon.

### How we measure it

Six metrics. The first is the one that matters for this project.

| Metric | What it asks | Ground truth | Direction | Source |
|---|---|---|---|---|
| **pFCS** | is the gene–gene **correlation structure** preserved? | answer-key cells | higher better | benchmark imputation task |
| MSE | how far is the guess from the real answer? | answer-key cells | lower better | benchmark imputation task |
| pDES | are differential-expression statistics preserved? | answer-key cells + labels | higher better | benchmark imputation task |
| gene-wise ρ | per-gene Pearson / Spearman | answer-key cells | 0.155–0.22 published ceiling | Polarbear / BABEL |
| densification | nonzero fraction of guessed rows vs real rows | real-row distribution | should not inflate | own script |
| downstream | 5 classification + 5 clustering metrics on the imputed data | labels | as the benchmark defines | benchmark imputation task |

**pFCS is the metric to lead with.** It is computed by taking the pairwise correlations of
the top 100 highly variable genes in the imputed data and in the held-out truth, then taking
the Pearson correlation of those two correlation matrices. It asks precisely the question
step 3 depends on: after imputation, do genes still correlate the way they really correlate?
MSE alone is insufficient, because a blurred matrix can have low MSE and still destroy
correlation structure — and correlation structure is the only thing the network model reads.

Definitions for MSE, pFCS and pDES follow the imputation-task evaluation in Liu et al.,
*Nat Methods* 2025 (10.1038/s41592-025-02856-3), which inherits pFCS from the Matilda
framework. Ready-made platforms: scIMC (Dai et al., *NAR* 2022, 10.1093/nar/gkac317) for imputation
benchmarking; the scTranslation benchmark (*Proc. ACM SIGKDD* 2026,
10.1145/3770855.3817464) for modality-translation metrics specifically.

**Negative control, required.** Repeat the whole measurement with randomized neighbours
(permute the reference). If the control scores nearly as well as the real imputation, the
metric is describing the data's general structure, not the imputation.

### Gate

Pass if all four hold:

1. pFCS clearly above the shuffled-reference control.
2. MSE clearly below the same control.
3. gene-wise ρ inside the published range — not near zero, and not suspiciously near 1.
4. The two blocks are on one scale: guessed rows are not markedly denser than measured
   rows, and per-gene mean and variance do not shift systematically between the two.

Failure signatures visible in the network, with no ground truth needed:

| Symptom | Cause |
|---|---|
| top edges are ribosomal / mitochondrial / HLA | correlation inflation (consequence 2) |
| edges present in both directions | symmetric smoothing |
| very large network with a flat ranking | same |

On failure: replace the imputer, or drop step 2 entirely and give step 3 the measured RNA
plus gene activity. That second option is legitimate — see the next section.

---

## 5. Step 3 — GRN inference

### What it does

A GRN is a table of TF → target links with a score. GRNBoost2 fits one regression per gene:
"which transcription factors predict this gene's expression?" It fits on whatever matrix it
is given, and it does not know which rows were guessed.

That last point is the reason the checks live upstream. A gradient-boosted tree finds any
pattern that exists, including the pattern "these 20 cells are a smoothed average", and it
reports no error while doing so.

### The files, and which of them is truth

Only the last two rows are ground truth, and they must never reach the inference step.

| File | What it is | Used for |
|---|---|---|
| `data/tf_only.txt` (829 symbols) | human transcription factors — the candidate **regulators**. Only these can appear on the left of an edge. | step 3 input |
| `data/trrust_tf.txt` (2,862 symbols) | every gene that acts as a regulator somewhere in TRRUST. A larger candidate set, and **not** a clean TF list: most of its entries are not transcription factors. | optional input |
| the matrix columns | the candidate **targets**, decided by the target rule below | step 3 input |
| `data/ground_truth/PBMC-TRRUST.csv` (9,720 `TF,TARGET` rows, no score) | a binary edge set for these cells | **evaluation only** |
| `data/ground_truth/PBMC-Blood.csv` (98,338 `TF,TARGET` rows) | the same, from whole blood | **evaluation only** |

The rule that follows from this table: nothing from `data/ground_truth/` may be read by the
inference step. It is a truth set, not a prior. If it leaks in — even indirectly, for example
by choosing the target genes to favour known edges — the evaluation is circular.

### Input

| Item | Form | Why it is needed |
|---|---|---|
| `X.npy` + `mask.npy` | cells × genes | from step 2 |
| regulator list | curated TF list, count recorded | what may appear on the left of an edge |
| target rule | stated rule, count recorded | what may appear on the right |
| cell grouping | one network per cell type, or pseudobulk | a single network over mixed cell types is dominated by differences between cell types |

Both counts are written down. Two runs with different target counts produce edge sets that
cannot be compared, and the difference is invisible unless the counts are recorded.

On the grouping choice: "Population-level comparisons of gene regulatory networks modeled on
high-throughput single-cell transcriptomics data", *Nat Comput Sci* 2024
(10.1038/s43588-024-00597-5), argues that single-cell-level networks are poorly
reproducible and that population-level comparison is the defensible unit.

### Five phrases in step 3

| Phrase | What it actually is |
|---|---|
| **regulators** ("TF list") | the genes allowed on the **left** of an edge. In this repo that is `tf_only.txt`. A regulator that is not a column of the matrix is dropped silently, so the list and the matrix have to be checked against each other once. |
| **target rule** | the rule that decides which genes are **columns of the matrix handed to GRNBoost2** — and therefore which genes get their own regression fitted. It matters because the public `grnboost2()` call takes no target list: it fits one regression per column, so the column set *is* the target list. A rule looks like "genes present in the matrix that are either a top-2000 highly variable gene or a TRRUST regulator", and its count is recorded. It is a scope-and-cost decision, not a claim about biology, and changing it changes every downstream edge count. |
| **cells grouped by cell type** | the **rows** of `X`, split by a barcode → cell-type table. One network is then inferred per group. The alternative is pseudobulk: average the cells within each cell type (optionally per donor) into a few hundred pseudo-cells and infer on those. |
| **importance** (the score) | how much predictive gain that TF contributed across all the trees that predict that target. It is a **ranking** score: not a p-value, not a correlation, and unsigned — it says nothing about activation versus repression. Compare importances within one run, not between runs or methods. |
| **three seeds** | GRNBoost2 is stochastic: each split considers only 10% of the regulators (`max_features=0.1`) and each tree sees 90% of the cells (`subsample=0.9`). Two runs on identical input therefore give different networks. Three seeds means running the same pipeline three times and measuring the overlap of the top-10k edges. That overlap is the noise floor; anything smaller than it is not a result. |

### Output

| Item | Form | Notes |
|---|---|---|
| `edges.csv` | TF, target, importance | one file per cell type |
| `run_info.txt` | cells, columns, seeds, workers, software versions | required for comparison |
| three seeds | the same run repeated | so seed-to-seed noise is measurable |
| top-N cap | e.g. top 50,000 edges by weight | the standard evaluation range |

The network is stored as an **edge list** — one row per link. That is the network, just not in
matrix form. The equivalent matrix is regulators × genes (2,827 × 4,389 ≈ 12.4 million cells
in this project), almost entirely zeros, so the edge list is both smaller and the form every
metric needs: AUPRC and early precision both take a ranked list of pairs. Pivot to a matrix
only for plotting.

Three seeds are not ceremony. The multimodal-GRN benchmark (GRETA) found these networks
highly sensitive to random seeds and to swapping a single pipeline step. If seed-to-seed
overlap is low, a difference between two methods is noise, and there is no way to know
without measuring it.

The 50,000 cap follows current benchmarking practice: geneRNIB evaluates only the top-ranked
TF–gene pairs, capped at 50,000 edges. A dense network's tail carries no information and
should not be scored.

### Methods

| Method | What it does | Choose it when | Reference |
|---|---|---|---|
| GRNBoost2 | boosting, one regression per gene | baseline | Moerman et al. 2019 |
| LINGER | paired RNA+ATAC, motif prior, atlas-scale bulk data as regularization | best accuracy is the goal; chromatin should stay in the model | Yuan et al., *Nat Biotechnol* 2025, 10.1038/s41587-024-02182-7 |
| SCENIC+ | enhancer→gene links from multiome, motif-based TF linking; outputs regulons | a stable evaluation unit is wanted | Bravo Gonzalez-Blas et al., *Nat Methods* 2023, 10.1038/s41592-023-01938-4 |
| HuMMuS | multilayer network over RNA, ATAC and protein interaction data; no imputation | mechanism-level view | *Bioinformatics* 2024, 10.1093/bioinformatics/btae143 |
| scMEGA | Bioconductor; peak–gene links then network | fast and inspectable, R-native | *Bioinformatics Advances* 2023, 10.1093/bioadv/vbad003 |

**The strategic point.** LINGER, SCENIC+ and HuMMuS consume chromatin directly. If one of
them works here, step 2 stops being a requirement and becomes an experiment that has to
justify itself. LINGER reports a four- to seven-fold relative accuracy increase over earlier
methods by combining paired multiome with a TF-motif prior and atlas-scale bulk data.

Keep the expectation honest: geneRNIB (below) finds that simple models with fewer assumptions
match or beat complex ones on this task. The boosting model is not the weak link.

### How we measure it

| Metric | What it asks | Ground truth | Tool |
|---|---|---|---|
| AUPRC | of all possible TF→target pairs, how well are the true ones ranked? | ChIP-seq / curated / knockout | BEELINE |
| early precision | precision among the top-ranked edges | as above | BEELINE |
| stability | Jaccard of top-10k edges across 3 seeds and subsamples | none needed | own script |
| topology | are hubs and in-degree structure recovered? | simulated networks | STREAMLINE |
| perturbation | does the network predict a measured knockout effect? | Perturb-seq | geneRNIB / OpenProblems |
| regulon activity | AUCell score of top regulons per cell type, and cell-type specificity | labels | SCENIC+/pySCENIC |

Two practical warnings. Use AUPRC rather than AUROC: with a few thousand true edges among
hundreds of thousands of possible pairs, AUROC looks flattering while AUPRC tells the truth.
And score only the top-ranked edges — see the 50,000 cap above.

Three truth tiers, following BEELINE (Pratapa et al., *Nat Methods* 2020,
10.1038/s41592-019-0690-6): synthetic networks with known structure (complete, artificial),
literature-curated Boolean models (realistic, small), and experimental data with ChIP or
knockout evidence (real, incomplete, noisy). Report all three. A method that wins only on
synthetic data has not won. STREAMLINE (Stock et al., *Bioinformatics* 2024,
10.1093/bioinformatics/btae267) adds the topological view — whether hubs and degree
distributions are recovered, which matters for robustness questions that edge-level AUPRC
does not capture.

geneRNIB ("a living benchmark for gene regulatory network inference", bioRxiv 2025,
10.1101/2025.02.25.640181) contributes two methodological points worth copying: each metric
is validated on a dataset by comparing the best method's score against a **negative control**
and by requiring a coefficient of variation of at least 0.2 across methods — a metric that
cannot discriminate between methods is dropped rather than reported. GRETA
(Badia-i-Mompel et al., bioRxiv 2025, 10.1101/2024.12.20.629764) is the modular framework
that decomposes multimodal GRN inference into swappable steps behind a fixed harness; it found
that predictive evaluation succeeds where causal recovery from perturbation assays does not.

### Gate

Pass if all three hold:

1. Seed-to-seed stability is high enough that a difference between two methods exceeds it.
2. Lineage regulators appear in the top edges of the matching cell type: PAX5/EBF1 →
   MS4A1/CD79A in B cells; SPI1/CEBPA → LYZ/S100A8 in myeloid; TBX21/GNLY in NK/T;
   CD3D/CD3E in T cells.
3. AUPRC and early precision beat the negative control on cell-type-matched truth.

Stability comes first, because a difference smaller than the method's own seed noise is not
a measurement. On failure: do not tune the booster — reconsider the matrix, or switch to a
multiome-native method.

---

## 6. Step 4 — End-to-end validation on simulated data

The three gates above grade each step in isolation. This section grades the pipeline as a
whole, against a network that is known to be there.

scMultiSim (*Nat Methods* 2025, 10.1038/s41592-025-02651-0) generates RNA, ATAC, RNA velocity
and spatial data from a gene regulatory network that the user specifies, while modelling cell
identity, cell–cell interactions, chromatin accessibility and technical noise. It demonstrates
benchmarking of integration, GRN inference and cell–cell interaction inference.

The test:

| Step | Metric | Truth |
|---|---|---|
| 1 | FOSCTTM | the simulated pairing |
| 2 | MSE, pFCS | the expression that was hidden |
| 3 | AUPRC, early precision, hub recovery | the network that was planted |

If the pipeline cannot recover a network it was given, with known structure and controlled
noise, then a poor result on real data says nothing about the biology. If it can, the
numbers from simulation become the calibrated expectation for real data — which is exactly
what is missing when reading published performance figures.

For plausibility on a finished network rather than pipeline correctness: CellOracle
(Kamimoto et al., *Nature* 2023, 10.1038/s41586-022-05688-9) performs in-silico gene
perturbation and compares the predicted shift in cell identity with expectation.

---

## 7. How the steps fit together

| Gate | Protects | From what |
|---|---|---|
| 1 — FOSCTTM + scIB | step 2 | a lookup space in which "nearest neighbour" means nothing, so the average is meaningless |
| 2 — pFCS + MSE + ρ | step 3 | a matrix that is confidently wrong: guessed rows are structured, so the network is structured and wrong |
| 3 — stability + lineage + AUPRC | step 4 | validating noise |
| end-to-end — simulation | the design itself | building four steps that were never shown to work together |

Four scripts in a row is a sequence. A pipeline is four steps where each output has been
shown good enough to be the next step's input. The gates are the difference, and each one
prevents one specific, predictable failure.

---

## 8. Run order and cost

| Order | Work | Cost |
|---|---|---|
| 1 | simulated data, all three grades | days, CPU |
| 2 | step 1 on real data, gate 1 | hours |
| 3 | step 2 on real data, gate 2 | hours to a day |
| 4 | compare two or three imputers, choose one | hours to days |
| 5 | step 3, once: three seeds, per cell type | the expensive part |
| 6 | step 4 on cell-type-matched truth | days |

Everything except step 3 runs on CPU in hours. Step 3 runs once, after its input has passed
a gate — not as the instrument used to test the earlier steps. That is the practical
consequence of the gating rule: the slow step stops being the debugging tool.

---

## 9. What good looks like

Set before running anything, so the result cannot be reinterpreted afterwards.

| Metric | A good result |
|---|---|
| FOSCTTM | 0.14–0.22 (published range); 0.5 is random; 0 is perfect |
| gene-wise ρ, ATAC→RNA | 0.155–0.22; not 0.9 — that suggests leakage |
| pFCS | high, and clearly above the shuffled-reference control |
| MSE | clearly below the same control |
| densification | guessed rows no denser than real rows |
| network stability | reported differences exceed seed-to-seed noise |
| AUPRC / early precision | above the negative control, on cell-type-matched truth |
| lineage regulators | present in the top edges of every cell-type network |

---

## 10. Sources

| Step | Reference | Identifier |
|---|---|---|
| 1 | Luecken et al., scIB, *Nat Methods* 2022 | 10.1038/s41592-021-01336-8 |
| 1 | Liu et al., multitask benchmarking, *Nat Methods* 2025 | 10.1038/s41592-025-02856-3 |
| 1 | Benchmarking single-cell multi-modal data integrations, *Nat Methods* 2025 | 10.1038/s41592-025-02737-9 |
| 1 | Qi et al., scGT, *Bioinformatics* 2025 | 10.1093/bioinformatics/btaf357 |
| 1 | scMODAL, *Nat Commun* 2025 | 10.1038/s41467-025-60333-z |
| 1 | SuperMap, *PNAS* 2026 | 10.1073/pnas.2505182123 |
| 1 | scHPGT, *Bioinformatics* 2026 | 10.1093/bioinformatics/btag696 |
| 1–2 | Ashuach, Gabitto et al., MultiVI, *Nat Methods* 2023 | 10.1038/s41592-023-01909-9 |
| 1 | FOSCTTM: definition and reference values (MMD-MA / UnionCom, SNARE-seq) | PMC8095090 |
| 2 | scGALA, *Nat Commun* 2025 | 10.1038/s41467-025-66644-5 |
| 2 | scDCT, *Brief Bioinform* 2025 | 10.1093/bib/bbaf400 |
| 2 | Cisformer, *Genome Biol* 2025 | 10.1186/s13059-025-03823-z |
| 2 | scButterfly, *Nat Commun* 2024 | 10.1038/s41467-024-47418-x |
| 2 | Zhang et al., Polarbear, *J Comput Biol* 2022 (also the ATAC→RNA Pearson 0.22 vs BABEL 0.155 figure) | 10.1089/cmb.2022.0264 |
| 2 | Dijkstra et al., imputation and network reconstruction, *Patterns* 2022 | 10.1016/j.patter.2021.100414 |
| 2 | scIMC, imputation benchmarking platform, *NAR* 2022 | 10.1093/nar/gkac317 |
| 2 | scTranslation benchmark, *Proc. ACM SIGKDD* 2026 | 10.1145/3770855.3817464 |
| 3 | Pratapa et al., BEELINE, *Nat Methods* 2020 | 10.1038/s41592-019-0690-6 |
| 3 | Stock et al., STREAMLINE, topological benchmarking, *Bioinformatics* 2024 | 10.1093/bioinformatics/btae267 |
| 3 | Yuan et al., LINGER, *Nat Biotechnol* 2025 | 10.1038/s41587-024-02182-7 |
| 3 | Bravo Gonzalez-Blas et al., SCENIC+, *Nat Methods* 2023 | 10.1038/s41592-023-01938-4 |
| 3 | HuMMuS, *Bioinformatics* 2024 | 10.1093/bioinformatics/btae143 |
| 3 | scMEGA, *Bioinformatics Advances* 2023 | 10.1093/bioadv/vbad003 |
| 3 | Population-level comparisons of GRNs, *Nat Comput Sci* 2024 | 10.1038/s43588-024-00597-5 |
| 3–4 | Badia-i-Mompel et al., GRETA, bioRxiv 2025 | 10.1101/2024.12.20.629764 |
| 3–4 | Nourisa et al., geneRNIB, bioRxiv 2025 | 10.1101/2025.02.25.640181 |
| 3–4 | Large-scale benchmark, single-cell perturbation, *Commun Biol* 2025 | 10.1038/s42003-025-07764-y |
| 4 | scMultiSim, *Nat Methods* 2025 | 10.1038/s41592-025-02651-0 |
| 4 | Kamimoto et al., CellOracle, *Nature* 2023 | 10.1038/s41586-022-05688-9 |
| 4 | Peidli et al., scPerturb, *Nat Methods* 2024 | 10.1038/s41592-023-02144-y |
