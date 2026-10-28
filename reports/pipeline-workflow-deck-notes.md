# Deck notes: from multiome data to a gene regulatory network

Companion to `reports/pipeline-workflow-design.pptx` (34 slides) and to
`reports/pipeline-workflow-design.md`.

Each slide below gives the text that goes **on the slide**, followed by the **speaker notes**.
Citations are kept on the slide, in short form; the full DOIs are in section 10 of the design
document.

Deck order: title → vocabulary → workflow at a glance → for each step (input & output,
methods, how good is it) → validation of the pipeline itself → gates → run order → what good
looks like → summary.

---

## Slide 1 — From multiome data to a gene regulatory network

*PIPELINE DESIGN*

**On the slide**

Integration  →  Imputation  →  GRN inference  →  Validation

For every step:  what goes in  ·  what comes out  ·  which methods  ·  how we judge the result

[ name ]     [ lab ]     [ date ]

**Speaker notes**

> Design talk, not a results talk. For each step: what goes in, what comes out, which methods, and how we judge the result.

---

## Slide 2 — Ten words, used throughout

**On the slide**

Nothing else in this deck is assumed

| Term | Plain meaning |
|---|---|
| RNA (scRNA-seq) | which genes each cell is transcribing |
| ATAC (scATAC-seq) | which DNA is open — what could be transcribed |
| multiome | RNA and ATAC measured in the same cell |
| modality | one kind of measurement (RNA is one, ATAC is another) |
| embedding | about 30 numbers per cell, so cells can be compared by distance |
| integration | one shared embedding for cells measured different ways |
| imputation | filling in a measurement a cell does not have |
| gene activity | a per-gene score computed from open chromatin near that gene |
| GRN | a table of TF → target gene links, with a score |
| regulon | a TF plus the genes it drives, scored per cell |

**Speaker notes**

> Ten words. Everything after this slide uses only these.

---

## Slide 3 — The workflow at a glance

**On the slide**

Each step's output is literally the next step's input

**1  INTEGRATION**

**IN**
- per-dataset PCA (50 PC)
- modality label
- dataset label
- paired cells set aside (answer key)

**OUT**
- H: cells × 30 numbers
- row-range map per dataset
- no genes, no imputation

**2  IMPUTATION**

**IN**
- H  (from step 1)
- measured RNA, one gene axis
- which cells need values
- answer-key cells excluded

**OUT**
- X: cells × genes
- mask: measured / guessed
- one normalization

**3  GRN INFERENCE**

**IN**
- X  (from step 2)
- TF list (regulators)
- target-gene rule
- cells grouped by type

**OUT**
- TF, target, importance
- cells used, 3 seeds
- one network per cell type

**4  VALIDATION**

**IN**
- edge table (step 3)
- independent truth
- a negative control

**OUT**
- AUPRC / early precision
- stability across seeds
- perturbation agreement
- a verdict

Only two objects matter — H and X — and both can be inspected without running a network.

**Speaker notes**

> Read it as a factory line: each output is the next step's input. Only two objects matter, H and X, and both can be inspected without running a network. Three terms in step 2 need defining once. Measured RNA is the gene-expression matrix of the cells that actually have RNA - the multiome cells - and it is the only source of the values copied onto ATAC cells. Which cells need values is the list of row numbers belonging to ATAC-only cells, because those are the rows being filled in. And the answer-key cells are a set of cells measured both ways, whose RNA we hide from the imputer and compare against afterwards; frozen means fixed once and never used to fit anything, so the comparison is honest.

---

## Slide 4 — Step 1 — Integration: input and output

**On the slide**

Two measurements of the same cells must land on one shared map

|  | Item | Shape / form | Why it matters |
|---|---|---|---|
| IN | per-dataset PCA | cells × 50, one file per dataset | integration never sees genes or peaks |
| IN | modality label | RNA / ATAC | the method must know which side a cell came from |
| IN | dataset label | one string per dataset | required by the batch-mixing metrics |
| IN | paired cells, held out | cell index, excluded from fitting | the only ground truth for this step |
| OUT | H.npy | cells × 30, float32 | the only thing step 2 consumes |
| OUT | blocks.json | dataset → [start, end) | row order must never be re-derived |
| OUT | per-dataset slices | one file each | saved once |
| NEVER | genes, peaks, imputation, reference choice | — | putting these here makes every later change unattributable |

**NOTE** — the ATAC 'feature space' is whether ATAC enters as peaks or as gene activity — it decides which integrators can run at all; fix it before comparing, and never compare scores across the two

**Speaker notes**

> Integration never sees genes or peaks, only 50 summaries per cell. The paired cells are held out from the start, because they are how this step gets graded. One input needs explaining: the ATAC feature space. ATAC is a matrix of genomic intervals - peaks - and different datasets almost never share the same intervals. Peaks keep full information but share no feature with RNA, so the integrator has to work without a shared axis, which is what the transport and manifold-alignment methods do. The alternative is gene activity: sum the open chromatin over each gene's body and promoter to get one number per gene, after which ATAC has a gene axis and any shared-axis method will run. It is lossy, because distal enhancers outside the window are dropped, and it correlates with expression, which flatters the scores. So record which one was used, and do not compare FOSCTTM across the two.

---

## Slide 5 — Step 1 — integration methods: what each one gives you

**On the slide**

Aligned embedding, latent, joint latent, shared factors — different words for the same object

| Method | What you get out | Citation |
|---|---|---|
| scSAGA | joint embedding (cells × 30) + the matching it was computed from | Aluru lab; bioRxiv 2026 |
| Pamona | joint embedding + matching; also labels cells shared vs dataset-specific | Cao et al., Bioinformatics 2021 |
| SCOT | matching + one dataset rewritten in the other's feature space; no low-dimensional embedding | Demetci et al., RECOMB 2020 |
| LIGER | joint embedding (cells × k) + per-dataset gene loadings (genes × k); ships its own imputeKNN | Welch et al., Cell 2019 |
| Seurat v5 | joint embedding; can also produce predicted expression for the other modality | Hao et al., Cell 2021 |
| Harmony | joint embedding; nothing else | Korsunsky et al., Nat Methods 2019 |

EVERY OUTPUT IS ONE OF THREE  a JOINT EMBEDDING (k numbers per cell — what step 2 needs)  ·  a MATCHING (cell-to-cell, convert it first)  ·  an IMPUTED MATRIX (genes × cells — step 2 is already done)

LIGER, Seurat and Harmony need RNA and ATAC on the same feature axis, so ATAC must enter as gene activity. scSAGA, Pamona and SCOT do not.

**Speaker notes**

> The point of this table is the middle column. These methods do not all return the same kind of object, and the papers use different words for it - aligned embedding, latent, joint latent, shared factors, corrected reduction. They all mean the same thing: k numbers per cell, every dataset in the same space, so distances between cells are meaningful. We call that a joint embedding. Everything reduces to three outputs: a joint embedding, which step 2 uses directly; a matching, which is a cell-to-cell probability matrix and has to be turned into a neighbour list or an embedding first; or an imputed expression matrix, in which case step 2 is already done. scSAGA gives a joint embedding of 30 dimensions plus the matching it came from. Pamona gives both, and also labels cells as shared or dataset-specific. SCOT is the one that is easy to get wrong: it gives a matching, plus a barycentric projection that rewrites one dataset in the other dataset's feature space - not a low-dimensional embedding. That matters twice. It is not a drop-in for a kNN imputer, because there is no latent space to measure distance in; and if the target domain is RNA, the projection is already an expression matrix, so SCOT can feed step two as a starting matrix or step three directly, as long as you say which. LIGER gives a joint embedding plus per-dataset gene loadings. Seurat gives a joint embedding and can also predict the other modality. Harmony gives an embedding and nothing else. One requirement is not in the table: LIGER, Seurat and Harmony need RNA and ATAC on the same feature axis, so ATAC must arrive as gene activity, whereas the transport methods do not care.

---

## Slide 6 — Step 1 — integration methods: multiome-native and recent

**On the slide**

Same three kinds of output; these differ in what prior they use

| Method | What you get out | Citation |
|---|---|---|
| MultiVI | joint embedding + an imputed expression matrix for the missing modality | Nat Methods 2023 |
| GLUE / scGLUE | joint embedding, built with a peak–gene prior; can also translate one modality into the other | Cao & Gao, Nat Biotechnol 2022 |
| scGT | joint embedding | Bioinformatics 2025 |
| scMODAL | joint embedding, built from a few known linked features | Nat Commun 2025 |
| SuperMap | a mapping between modalities learned from unpaired data; no joint embedding by default | PNAS 2026 |
| scHPGT | joint embedding, built with a regulatory prior | Bioinformatics 2026 |

**READ FIRST** — Liu et al., Nat Methods 2025 — 40 methods benchmarked on dimension reduction, batch correction, clustering, label transfer and imputation, with a decision guideline

**Speaker notes**

> Two kinds of method are missing from the previous table. MultiVI is the one that does steps one and two in a single model, so the joint embedding and the imputed expression come from the same fit. GLUE aligns the modalities through a prior of peak-to-gene links, which makes it prior-driven rather than anchor-free. Everything else here is the 2025-26 generation, and in every case what you get out is a joint embedding. Before choosing from either table, read the benchmark strip: Liu et al. 2025 evaluates forty methods on exactly these tasks and publishes a decision rule.

---

## Slide 7 — Step 1 — Integration: how good is it?

**On the slide**

Four metrics; the first needs paired cells and is the strictest

| Metric | What it asks | Truth needed | Good looks like |
|---|---|---|---|
| FOSCTTM | for each cell, what fraction of other cells are closer to it than its true cross-modal partner is? | paired cells, held out | 0 = perfect, 0.5 ≈ random. Published: 0.14–0.16 (MMD-MA), 0.22 (UnionCom), 0.39 (unnormalized) |
| label conservation | are cell types still separated? ARI, NMI, cLISI, silhouette, isolated-label F1, graph connectivity | cell-type labels | close to unintegrated labels, far from a shuffled control |
| batch correction | are datasets mixed? kBET, iLISI, batch silhouette | dataset IDs | improves while the label metrics hold |
| label transfer | can a label learned on RNA be assigned to the ATAC half? | paired cells + labels | > 0.85 on major lineages |

**PASS IF** — FOSCTTM below the shuffled-partner control  +  label transfer > 0.85  +  batch metrics up while label metrics hold.   Not a gate: scSAGA's 'global alignment score' — it moves with the anchor.

scib implements the middle two rows (14 metrics). Report both halves separately — a single combined score hides the batch-versus-biology trade (Luecken et al., Nat Methods 2022).

**Speaker notes**

> FOSCTTM is simple: for a cell with a known partner, what fraction of other cells are closer to it than that partner is. Perfect is 0, random is about 0.5, and the best published numbers on real paired data are 0.14 to 0.22 - so the scale is calibrated before we use it. Report the scIB label and batch halves separately; a single combined score can hide a method that removed batch by destroying biology.

---

## Slide 8 — Step 2 — Imputation: input and output

**On the slide**

One gene axis, one normalization, one frozen held-out set

|  | Item | Shape / form | Why it matters |
|---|---|---|---|
| IN | H.npy + blocks.json | cells × 30 | the lookup space from step 1 |
| IN | measured RNA | cells × genes, ONE gene axis | the values being copied onto ATAC cells |
| IN | normalization | declared, e.g. log1p CPM/1e4 | measured and guessed rows must share a scale |
| IN | which cells need values | row indices of ATAC-only cells | the rows of X being filled in |
| IN | answer-key cells | paired cells, fixed and excluded | hidden RNA we compare against; if they leak in, the test is void |
| OUT | X.npy | cells × genes, float32 | the only input to step 3 |
| OUT | mask.npy | cells × genes, bool | per cell per gene: measured or guessed |
| OUT | genes.txt, norm.txt | one axis, one normalization | written once, never re-derived |
| OUT | reference set used | cell indices | proof the answer key was never inside |

**WHY THE MASK** — without the mask nobody downstream can point at a row and say 'that one is an average' — so nobody can ask whether step 2 helped

**Speaker notes**

> Three constraints: one gene axis, one normalization, one frozen answer key. Two terms worth spelling out. The mask is a true/false grid the same shape as X, true where the value was measured and false where it was guessed; it is how any later step can ask whether a number is real, and it is what lets step 3 run with and without the guessed rows. One normalization means exactly one declared transform, for example log1p of counts per 10k written to norm.txt, applied so that measured and guessed rows sit on the same scale - if the imputer outputs something else, such as a probability between 0 and 1, it has to be mapped into that space. Otherwise the two blocks are separable and the network learns which block a cell came from instead of biology. The quick check is per gene: compare mean, variance and nonzero fraction between measured and guessed rows. They should be close; a systematic shift means the scales do not match.

---

## Slide 9 — Step 2 — imputation from H: what actually exists

**On the slide**

The weights are available off the shelf — the transfer is not

| Tool | What it consumes | What you actually get | Citation |
|---|---|---|---|
| reverse-imputeKNN  (this repo) | H + the measured RNA + a reference block | ATAC→RNA kNN average — the transfer is implemented | run.py |
| LIGER  imputeKNN | H (its H.norm) + the reference's peaks | RNA→PEAKS only: modal = "atac" is hard-coded, so the direction is wrong for us | Welch et al., Cell 2019 |
| moscot | a joint latent passed as joint_attr | a transport plan — the transfer is your code | Klein et al., Nature 2025 |
| SCOT | no latent needed | a matching + projection into feature space — no latent space | Demetci et al., RECOMB 2020 |
| Symphony  mapQuery | reference embedding + cluster centroids | query cells in the reference space + labels — transfer is yours | Kang et al., Nat Commun 2021 |
| MOFA+  project_data | factors + gene loadings | factors for new cells — reconstruction is yours | Argelaguet et al., Genome Biol 2020 |

THE GAP  nothing published takes H and returns gene expression for ATAC cells: LIGER's imputer runs the opposite direction, and the OT tools give you the weights, not the average. That last mile is one line of code — and it is what reverse-imputeKNN is

**Speaker notes**

> This is the honest inventory, and it is the most useful slide in the section because it answers what we can actually run. Nothing published takes H and hands back gene expression for ATAC-only cells. LIGER is the closest and it is the wrong direction: imputeKNN is hard-wired to impute PEAKS for query cells - it pulls rawPeak from the reference, forces the query dataset to modal equals atac, and then normalises peaks. Gene to peak, not peak to gene. That is exactly why our step is called reverse imputeKNN: it is the transpose of LIGER's function, which did not exist. The other tools give you the weights but not the transfer. moscot returns a transport plan in a shared latent. SCOT returns a matching and a projection in feature space, with no latent at all. Symphony returns query cells placed in the reference space plus labels. MOFA+ returns factors. In every case the weighted average over reference expression is code you write - and it is short. So the real design choice for step two is not which package to install, it is which weighting rule to use: kNN with exp minus distance, plain kNN, or an optimal-transport plan; and which reference set. That is a much smaller and more tractable design space than the method table made it look.

---

## Slide 10 — Step 2 — imputation methods that build their own embedding

**On the slide**

These replace step 1: the integration happens inside the imputer, so there is no H to hand over

| Method | What it consumes | Kind of model | Citation |
|---|---|---|---|
| scGALA | measured RNA + ATAC as gene activity + its own cell alignment | graph attention network | Nat Commun 2025 |
| scDCT | raw RNA and ATAC matrices; paired or unpaired | conditional diffusion | Brief Bioinform 2025 |
| Cisformer | ATAC only, 10,000 open regions per cell; pretrained | transformer | Genome Biol 2025 |
| scButterfly | raw RNA and ATAC matrices; paired or unpaired | two VAEs + translator | Nat Commun 2024 |
| MultiVI | raw RNA and ATAC counts, with a mask for which modality each cell has | joint VAE | Nat Methods 2023 |
| Polarbear | raw RNA and ATAC matrices + co-assay cells for the translator | VAE + linear map | J Comput Biol 2022 |

THEY REPLACE STEP 1  easier to adopt and better on benchmarks, but they merge steps 1 and 2 — the step-1 gate then applies to the model's internal embedding, and there is no H file in the middle

Ceiling for this task: ATAC→RNA gene-wise Pearson ≈ 0.22 (Polarbear), 0.155 (BABEL), on SNARE-seq. 0.9 suggests leakage; 0.05 means nothing was produced.

**Speaker notes**

> These are the recent, better-benchmarked options, and they all share one property: they take the raw matrices and learn their own joint representation internally. scGALA needs reference RNA, ATAC converted to gene activity, and its own anchor file; scDCT and scButterfly train VAEs on raw matrices; Cisformer takes ATAC open regions only; MultiVI models both modalities with a mask; Polarbear adds a translator layer between two VAEs. Adopting any of them does not remove the need to integrate - it moves the integration inside the imputer. The cost is that there is no H artifact to gate on, so the step-1 check has to be applied to the model's internal embedding instead, and the pipeline becomes three steps rather than four. Ceiling to keep in mind: about 0.22 gene-wise Pearson for ATAC to RNA.

---

## Slide 11 — Step 2 — Imputation: what averaging does to the matrix

*The cheap option fails in a predictable way — this is the mechanism*

**On the slide**

**1  Zeros disappear** — An average of 20 cells is rarely zero, so guessed cells look richer than measured ones.

**2  Correlations inflate** — Two genes present in the same neighbourhood become correlated, because we averaged them, not because a TF links them.

**3  The block is detectable** — Guessed rows have a different distribution from measured rows, so a model can separate the two.

**4  The GRN chases 1–3** — A model that ranks genes by mutual predictability ranks the inflated correlations first.

**The visible signature of consequence 4** — a network whose strongest links are housekeeping genes (ribosomal, mitochondrial, HLA), because those are high and correlated in every neighbourhood.

Not hypothetical: Dijkstra et al., *Patterns* 2022 (10.1016/j.patter.2021.100414) measured networks before and after imputation and found inflated gene–gene correlations, altered network structure, and a possible decrease in reconstruction performance.

**Speaker notes**

> This is the mechanism, and it is why the metrics on the next slides matter. Averaging is smoothing, and smoothing manufactures correlation. Whatever truly regulates what, the smoothed matrix insists these genes always appear together, and the boosting model believes it. Housekeeping genes are high and correlated in every neighbourhood, so they surface at the top of the ranking.

---

## Slide 12 — Step 2 — Imputation: how good is it?

**On the slide**

Lead with pFCS — it is the one metric that matches what step 3 reads

| Metric | What it asks | Truth needed | Good looks like |
|---|---|---|---|
| pFCS | is the gene–gene CORRELATION structure preserved? Pearson of the two correlation matrices, top-100 HVG | held-out cells | high, and clearly above the control |
| MSE | how far is the guess from the real answer? | held-out cells | low, and clearly below the control |
| pDES | are differential-expression statistics preserved? (limma per cell type) | held-out cells + labels | high |
| gene-wise ρ | per-gene Pearson / Spearman | held-out cells | 0.155–0.22 (published) |
| densification | nonzero fraction, guessed rows vs measured rows | measured-row distribution | no inflation |
| downstream | 5 classification + 5 clustering metrics on the imputed data | labels | as the benchmark defines |

**NEGATIVE CONTROL** — repeat the measurement with randomized neighbours. If the control scores as well, the metric describes the data, not the imputation.

Definitions follow the imputation task in Liu et al., Nat Methods 2025; scIMC (NAR 2022) and the scTranslation benchmark (ACM SIGKDD 2026) are ready-made harnesses.

**Speaker notes**

> pFCS asks the exact question step 3 depends on: after imputation, do genes still correlate the way they really correlate? A blurred matrix can have low MSE and still destroy correlation structure, and correlation structure is the only thing the network model reads.

---

## Slide 13 — Step 3 — GRN inference: input and output

**On the slide**

The model does not know which rows were guessed

|  | Item | Shape / form | Why it matters |
|---|---|---|---|
| REQUIRED | X.npy + mask.npy | cells × genes | the matrix from step 2, guessed rows flagged |
| REQUIRED | regulator list | curated TFs; count recorded | what may appear on the left of an edge |
| DECIDE | target rule | default: every column; or a stated rule | the column set IS the target list — changes cost 8× |
| DECIDE | cell grouping | default: one network; or per cell type | a mixed network is dominated by cell-type differences |
| OUT | edges.csv | TF, target, importance | one file per cell type |
| OUT | run_info.txt | cells, columns, seeds, software versions | required for any comparison |
| OUT | three seeds | the same run repeated | makes seed-to-seed noise measurable |
| OUT | top-50k cap | ranked edges used for scoring | the standard evaluation range |

WHY RECORD  the two DECIDE rows are not needed to run the step, but they change the cost and the denominator of every metric — write them down, or two runs cannot be compared

**Speaker notes**

> The booster does not know which rows were guessed. Only two of these four rows are actually required: the matrix and the regulator list. The other two are decisions with defaults - if you say nothing, every column becomes a target and you get one network over all cells - so they are not inputs you must go and find, but they change the cost and the metric denominator enough that leaving them unrecorded makes the run incomparable. Four things worth spelling out. The regulators are the genes allowed on the left of an edge - in this repo tf_only.txt, 829 human transcription factors; a regulator that is not a column of the matrix is dropped silently. The target rule decides which genes are columns of the matrix handed to the booster, and because grnboost2 takes no target list and fits one regression per column, the column set IS the target list - it is a scope and cost decision, and changing it changes every edge count. Cells grouped by cell type means splitting the rows by a barcode-to-cell-type table and inferring one network per group, or pseudobulking within a cell type for reproducibility. And the output is a network stored as an edge list, one row per link: TF is the regulator on the left, target is the gene being predicted, importance is how much predictive gain that TF contributed to that gene's model - a ranking score, not a p-value, not a correlation, and unsigned. To be clear, the network IS a matrix - regulators by genes, with the weight of each edge in the cell. What is written to disk is that same matrix in sparse form, one row per non-zero entry, and the two convert in one line either way. The sparse form is the default not because the matrix is big - it is about 47 megabytes here - but because a zero in a dense matrix means both no edge and zero weight, every metric wants a ranked list of pairs, and importances are scaled per target so reading across a row compares things that are not comparable.

---

## Slide 14 — Step 3 — GRN inference: methods

**On the slide**

Three of these consume chromatin directly, which makes step 2 optional

| Method | What it does | Choose it when | Citation |
|---|---|---|---|
| GRNBoost2 | boosting; one regression per gene | baseline | Moerman et al. 2019 |
| LINGER | paired RNA+ATAC, motif prior, atlas-scale bulk data | best accuracy; keep chromatin in the model | Nat Biotechnol 2025 |
| SCENIC+ | enhancer→gene links from multiome; outputs regulons | a stable evaluation unit is wanted | Nat Methods 2023 |
| HuMMuS | multilayer RNA + ATAC + protein network; no imputation | mechanism-level view | Bioinformatics 2024 |
| scMEGA | peak–gene links, then network; Bioconductor | fast and inspectable, R-native | Bioinform Adv 2023 |

**The strategic point** — LINGER, SCENIC+ and HuMMuS use the ATAC directly. If one of them works here, step 2 stops being a requirement and becomes an experiment that has to justify itself. LINGER reports a 4–7× relative accuracy gain over earlier methods.

Keep it honest: geneRNIB (2025) finds simple models with fewer assumptions match or beat complex ones here — the booster is not the weak link.

**Speaker notes**

> This is the decision worth discussing with the group. LINGER, SCENIC+ and HuMMuS take chromatin as chromatin. If one works here, imputation becomes an experiment rather than a requirement.

---

## Slide 15 — Step 3 — GRN inference: how good is it?

**On the slide**

AUPRC, not AUROC; and only the top-ranked edges

| Metric | What it asks | Truth needed | Good looks like |
|---|---|---|---|
| AUPRC | of all possible TF→target pairs, how well are the true ones ranked? | ChIP-seq / curated / knockout | above the negative control |
| early precision | precision among the top-ranked edges | as above | BEELINE's headline metric |
| stability | Jaccard of top-10k edges across 3 seeds and subsamples | none — measures our own noise | reported differences must exceed it |
| topology | are hubs and degree structure recovered? | simulated networks (STREAMLINE) | comparable to the planted network |
| perturbation | does the network predict a measured knockout effect? | Perturb-seq | positive agreement |
| regulon activity | AUCell score of top regulons per cell type | cell-type labels | cell-type specific, not global |

**PASS IF** — differences must exceed seed noise  +  lineage TFs in the top edges (PAX5/EBF1→MS4A1 in B; SPI1/CEBPA→LYZ in myeloid; TBX21/GNLY in NK)  +  AUPRC above control

Three truth tiers, following BEELINE: synthetic networks (complete, artificial), literature Boolean models (realistic, small), experimental ChIP / knockout (real, incomplete). Report all three. geneRNIB additionally gates each metric by its own negative control and a coefficient of variation ≥ 0.2 across methods.

**Speaker notes**

> Use AUPRC rather than AUROC: with a few thousand true edges among hundreds of thousands of possible pairs, AUROC looks flattering while AUPRC tells the truth. Score only the top-ranked edges, and report stability before anything else. On the three seeds: the booster is stochastic, since each split looks at ten percent of the regulators and each tree at ninety percent of the cells, so the same input gives a slightly different network every run. Three seeds, then overlap of the top ten thousand edges, gives the noise floor - anything smaller than that is not a result.

---

## Slide 16 — Six ways to get from multiome data to a network

**On the slide**

Only way 1 keeps a separate, gateable embedding — the others buy scale or justification

| Way | What it produces | What the GRN sees | Scales to | Main risk |
|---|---|---|---|---|
| 1  integrate, then impute from H | a joint embedding H, then a gene × cell matrix | measured + imputed RNA, all cells | OT to 1M cells | the last-mile transfer is bespoke code |
| 2  one joint model does both | imputed RNA straight from raw RNA + ATAC | measured + imputed RNA, all cells | minibatch training | no H to gate; ATAC→RNA ceiling ≈ 0.22 |
| 3  gene activity, no translation | ATAC projected onto the gene axis | measured RNA + gene activity, all cells | sparse matmul, linear | mixes two different quantities in one matrix |
| 4  multiome-native network | TF→target edges with enhancer support | peaks + motifs, no gene × cell matrix | ~10⁵ cells in practice | needs motifs + peak–gene links; heavy |
| 5  aggregate first | metacells or pseudobulk groups | few thousand pseudo-cells | collapses 1M → 10³–10⁴ | loses per-cell resolution |
| 6  TF activity as the regulator side | per-cell TF motif activity from ATAC | TF activity + measured RNA targets | sparse, linear | activity is a motif proxy, not expression |

**Speaker notes**

> This is the survey part of the talk, and the point is that the architecture is a choice, not a given. Six routes exist in the literature to get from multiome data to a network. Only the first produces a separate embedding you can inspect and grade. The others trade that for scale, for accuracy, for a better-justified network, or for avoiding imputed values altogether. The columns that matter are the third and fourth: what the GRN actually sees, and how far the route scales. At a million cells the fourth column is what decides, and most of these routes scale without any optimal transport at all.

---

## Slide 17 — Way 1 — integrate, then impute from H: what goes in and what comes out

**On the slide**

The architecture in this deck: H is an artifact, and both steps can be graded

**integrate** — OT / manifold alignment

**weight** — kNN or OT over reference

**impute** — reference RNA @ weights

**GRN** — one network, all cells

| Step | What goes in | What comes out |
|---|---|---|
| 1 integrate | per-dataset PCA, modality + dataset labels, paired cells held out | H (cells × 30) + block map; no genes |
| 2 weight | H + measured RNA + which rows need values | weights (query × reference cells), sparse |
| 3 impute | weights + measured RNA on one gene axis, one normalization | X (cells × genes) + mask + genes.txt |
| 4 GRN | X + mask, TF list, target rule, cells grouped by type | edges: TF, target, importance |

THE GRN RECEIVES  measured RNA rows plus weighted averages for the ATAC rows, with a mask saying which is which

**Speaker notes**

> Step by step, way 1 is: an embedding from per-dataset PCA, then weights from a neighbour or transport lookup, then a weighted average of measured RNA, then the network. Four artifacts, four chances to grade. The important line is the bottom strip: what the GRN receives is a matrix where some rows are measured and some are averages, and the mask is what keeps that distinction alive.

---

## Slide 18 — Way 1 — methods, and how the route is judged

**On the slide**

The alignment is off the shelf; the transfer is not

| Method | What it does | Citation |
|---|---|---|
| scSAGA | optimal transport with sampled anchors → H plus the plan | bioRxiv 2026 |
| Pamona | partial Gromov–Wasserstein coupling + a common space | Bioinformatics 2021 |
| SCOT | GW coupling + barycentric projection into feature space | RECOMB 2020 |
| kNN / OT weighting | the transfer itself — no published implementation for this direction | run.py |

| Metric | What it measures |
|---|---|
| FOSCTTM | is each cell nearest its true cross-modal partner? 0 = perfect, 0.5 = random |
| pFCS | is gene–gene correlation preserved after averaging? the metric step 3 reads |
| seed / k stability | does the network survive 3 seeds and k = 10 vs 30? |

WHAT TO WATCH  the transfer is untested code by construction — no benchmark and no reference implementation, so its pFCS gate against shuffled neighbours is the only evidence that it works

**Speaker notes**

> The alignment half is well served: scSAGA, Pamona and SCOT all give you an embedding and a plan. The transfer half is not served at all - that row is our own code. Which is why the validation column matters more here than anywhere else: FOSCTTM grades the embedding, pFCS grades the average, and seed stability grades the network. There is no published benchmark standing behind the middle row, so the shuffled-neighbour control is the only external check we get.

---

## Slide 19 — Way 2 — one joint model: what goes in and what comes out

**On the slide**

Integration happens inside the model, so there is nothing in the middle to hand over

**fit** — VAE / diffusion, raw counts

**latent** — the model's own space

**decode** — RNA for ATAC cells

**GRN** — one network, all cells

| Step | What goes in | What comes out |
|---|---|---|
| 1 fit | raw RNA counts + ATAC counts, plus a mask of which modality each cell has | a trained model; nothing else yet |
| 2 latent | the fitted model | z (cells × d) — internal; export it if you want a gate |
| 3 decode | z + the decoder | X for ATAC cells, stacked with measured RNA |
| 4 GRN | X + mask, TF list, target rule | edges: TF, target, importance |

THE GRN RECEIVES  the same kind of matrix as way 1, but produced in one fit rather than two

**Speaker notes**

> Same destination, one fewer artifact. You feed raw counts with a modality mask, get a latent, decode RNA for the cells that lacked it, and stack. Nothing in the middle is inspectable unless you deliberately export the latent, which is why the second row says so explicitly - if you do not save it, there is nothing to gate.

---

## Slide 20 — Way 2 — methods, and how the route is judged

**On the slide**

Best published accuracy for imputation; least attributable

| Method | What it does | Citation |
|---|---|---|
| MultiVI | joint VAE; masks the missing modality per cell | Nat Methods 2023 |
| GLUE / scGLUE | peak–gene guidance graph; adversarially aligned latents | Nat Biotechnol 2022 |
| scGALA | GAT mapper on gene activity + its own anchor pairs | Nat Commun 2025 |
| scDCT / scButterfly / Polarbear | diffusion, dual VAE, VAE + linear translator | 2025 / 2024 / 2022 |

| Metric | What it measures |
|---|---|
| held-out MSE, pFCS | mask paired cells, impute them, compare with the hidden RNA |
| FOSCTTM on the exported latent | possible only if z is saved; otherwise no step-1 gate |
| run-to-run agreement | two seeds give different latents — measure edge overlap |

WHAT TO WATCH  the latent is not identifiable: retraining reorders it, so the step-1 gate must be recomputed per fit and hyperparameters move the result more than the data does

**Speaker notes**

> These are the best-benchmarked imputers and they are all deep generative models: a joint VAE, a graph-guided VAE, a graph attention mapper, a diffusion model, a dual VAE, a linear translator on top of a VAE. The judgement is the same three checks as way 1 provided you export the latent, plus one that is specific to this route: two seeds give two different latent spaces, so reproducibility has to be measured on the edges, not on the metric.

---

## Slide 21 — Way 3 — gene activity: what goes in and what comes out

**On the slide**

No translation at all: ATAC is projected onto genes and stacked

**score** — ATAC → per-gene score

**stack** — RNA + gene activity

**mask** — flag scored rows

**GRN** — one network, all cells

| Step | What goes in | What comes out |
|---|---|---|
| 1 score | ATAC peaks + gene annotation (bodies, promoters) | gene activity (genes × cells) |
| 2 stack | gene activity + measured RNA on one shared gene axis | X (cells × genes) |
| 3 mask | row provenance | mask: measured vs scored |
| 4 GRN | X + mask, TF list, target rule | edges: TF, target, importance |

THE GRN RECEIVES  measured RNA plus accessibility-derived scores for the ATAC rows — a different quantity, on the same gene axis

**Speaker notes**

> This route deletes the hard step instead of solving it. ATAC goes to per-gene scores using the annotation, those scores are stacked with measured RNA on the same gene axis, the mask records which rows are scores, and the network runs on the whole matrix. Nothing is imputed, so nothing can be imputed wrongly - but the rows are not the same kind of number, and saying that out loud is the whole point of the mask here.

---

## Slide 22 — Way 3 — methods, and how the route is judged

**On the slide**

Cheapest route to scale, with the clearest semantic caveat

| Method | What it does | Citation |
|---|---|---|
| ArchR addGeneScoreMatrix | accessibility weighted across gene body + promoter | Nat Genet 2021 |
| Signac GeneActivity | promoter extended 2 kb; works from fragments or peaks | Nat Methods 2021 |
| + peak–gene links | add distal enhancers instead of a fixed window | ArchR addPeak2GeneLinks |

| Metric | What it measures |
|---|---|
| gene activity vs measured RNA | per-gene Pearson on paired cells — expect modest, not high |
| scale match | distribution of scored vs measured rows must not shift systematically |
| does it add or dilute? | GRN on mixed rows vs GRN on RNA-only — report the delta |

WHAT TO WATCH  a scored row is not a worse expression row, it is a different quantity: in one matrix the booster prefers whichever block has more variance

**Speaker notes**

> ArchR and Signac both compute gene scores; the difference is the window, and ArchR can add peak-gene links if you want distal enhancers rather than a fixed window. The judgement has to start with an honest measurement: correlate the scores against measured RNA on paired cells. Expect a modest correlation, not a high one. Then check that the two blocks have comparable distributions, and finally ask the only question that matters - does adding the scored rows help or dilute the network?

---

## Slide 23 — Way 4 — multiome-native network: what goes in and what comes out

**On the slide**

Chromatin stays chromatin; the network is built from peaks, motifs and links

**peaks** — ATAC + RNA

**links** — peak→gene, motif→TF

**network** — TF→target

**grade** — regulons

| Step | What goes in | What comes out |
|---|---|---|
| 1 peaks | paired RNA + ATAC, motif database, gene annotation | region accessibility + motif matches |
| 2 links | regions + annotation, optionally accessibility correlation | peak → gene links |
| 3 network | links + motifs + TF expression | TF→target edges with enhancer support |
| 4 grade | regulons + per-cell activity | eRegulons + AUCell scores |

THE GRN RECEIVES  TF→target edges with the enhancers that justify them — no gene × cell matrix at all

**Speaker notes**

> No gene by cell matrix anywhere. Peaks and motifs give region-level TF evidence, peak-gene linking turns that into enhancer support, and the output is a network plus regulons scored per cell. This is the only route where the chromatin evidence is inside the network rather than folded into a matrix, and it is the better artifact for a biology audience.

---

## Slide 24 — Way 4 — methods, and how the route is judged

**On the slide**

The only route where the ATAC evidence stays inside the network

| Method | What it does | Citation |
|---|---|---|
| SCENIC+ | topics → motif enrichment → enhancer links → eRegulons | Nat Methods 2023 |
| LINGER | paired multiome + motif prior + atlas-scale bulk | Nat Biotechnol 2025 |
| HuMMuS | multilayer RNA + ATAC + protein-interaction network | Bioinformatics 2024 |
| FigR / scMEGA | chromVAR activity + peak–gene correlation; R-native peak–gene route | Cell Genomics 2022 / Bioinform Adv 2023 |

| Metric | What it measures |
|---|---|
| ChIP / ReMap overlap | cell-type-matched AUPRC — not the whole-blood list |
| regulon specificity | AUCell per cell type; on-target vs off-target |
| perturbation | does the network predict a measured knockout? |

WHAT TO WATCH  peak–gene linking is correlation, not causality: the enhancer support is the weakest claim in the output and should be reported as a hypothesis

**Speaker notes**

> SCENIC+ is the reference implementation, LINGER and HuMMuS are the recent accuracy plays, FigR and scMEGA are lighter. Judgement is different in kind: you grade the network against cell-type-matched ChIP or ReMap, you grade the regulons by specificity across cell types, and if you are serious you grade causality against a perturbation. The honest caveat is the peak-gene link: it is a correlation, so the enhancer support is a hypothesis, not evidence.

---

## Slide 25 — Way 5 — aggregate first: what goes in and what comes out

**On the slide**

Change the unit of inference rather than the method

**aggregate** — metacells / pseudobulk

**10³–10⁴** — pseudo-cells

**infer** — any GRN method

**grade** — population truth

| Step | What goes in | What comes out |
|---|---|---|
| 1 aggregate | a cell × gene matrix from any route + a cell → group assignment | pseudo-cells (10³–10⁴ × genes) + labels |
| 2 infer | pseudo-cell matrix, TF list, target rule | edges: TF, target, importance |
| 3 grade | edges + cell-type truth | AUPRC / early precision |

THE GRN RECEIVES  a small matrix of pseudo-cells, which any GRN method handles — dropout is reduced and the network is far more reproducible

**Speaker notes**

> The smallest pipeline in the deck. Take any cell by gene matrix, plus an assignment of cells to groups, and aggregate. What comes out is a few thousand pseudo-cells, and everything downstream is cheap. This is the route that makes scale a non-issue, and it is also the route with the best evidence for reproducibility, because population-level networks are far more stable than single-cell ones.

---

## Slide 26 — Way 5 — methods, and how the route is judged

**On the slide**

The cheapest route to scale, and the least per-cell

| Method | What it does | Citation |
|---|---|---|
| SEACells | archetypal cells from a low-dimensional embedding, soft membership | Nat Biotechnol 2023 |
| MetaCell | balanced kNN graph partitions with a size target | Genome Biol 2019 |
| pseudobulk | sum counts per cell type × donor | standard practice |

| Metric | What it measures |
|---|---|
| metacell purity | each metacell should be one cell state — check the mixing |
| partition stability | re-run the aggregation; do you get the same units? |
| edge concordance | Jaccard against the single-cell network |

WHAT TO WATCH  the unit changes, so the numbers are not comparable to per-cell networks, and a rare population can disappear entirely

**Speaker notes**

> SEACells and MetaCell both build metacells, from a low-dimensional embedding and from a kNN graph respectively; pseudobulk is the simplest version and needs a donor column. The checks are about the aggregation itself: is each metacell one cell state, do you get the same units when you re-run, and how much does the final network agree with the single-cell one. The last one is the number that tells you what you gave up.

---

## Slide 27 — Way 6 — TF activity as the regulator side: what goes in and what comes out

**On the slide**

ATAC supplies the regulators, RNA supplies the targets

**activity** — chromVAR on ATAC

**TF × cell** — deviation scores

**link** — activity → target

**network** — TF→target

| Step | What goes in | What comes out |
|---|---|---|
| 1 activity | ATAC peaks + motif matches (chromVAR) | TF × cell deviation scores |
| 2 link | activity + measured RNA targets | TF→target regression |
| 3 network | coefficients or importances | edges: TF, target, weight |

THE GRN RECEIVES  a TF × cell activity matrix alongside measured expression — no gene values are imputed at any point

**Speaker notes**

> The asymmetric route: ATAC is used only for the regulator side, as per-cell motif activity, and RNA supplies the targets. There is no imputation at any point, which is the attraction if what you distrust is imputed expression. The output is a TF-by-cell activity matrix and a network regressed from activity onto targets.

---

## Slide 28 — Way 6 — methods, and how the route is judged

**On the slide**

Useful when imputed expression is what you distrust

| Method | What it does | Citation |
|---|---|---|
| chromVAR | motif accessibility deviations per cell | Nat Methods 2017 |
| FigR | chromVAR activity + peak–gene correlation → TF–gene network | Cell Genomics 2022 |
| SCENIC+ motif enrichment | an alternative per-cell TF activity estimate | Nat Methods 2023 |

| Metric | What it measures |
|---|---|
| activity vs TF expression | sanity: the same TF should agree between the two views |
| ChIP overlap | TF→target AUPRC, cell-type matched |
| knockout response | predicted vs measured perturbation |

WHAT TO WATCH  a TF with no motif is invisible, GC-rich motifs are over-called, and activity is accessibility rather than protein

**Speaker notes**

> chromVAR computes the activity, FigR is the pipeline that turns activity plus peak-gene correlation into a TF-gene network, and SCENIC+ motif enrichment is an alternative estimator. The first check is almost free and worth doing every time: does the activity of a TF agree with its own measured expression? Then ChIP overlap, then perturbation. The standing caveat is that a TF with no motif in the database is simply invisible to this route.

---

## Slide 29 — How to choose between the ways

**On the slide**

They compose — and the binding constraint decides

| If the binding constraint is… | Take | Because |
|---|---|---|
| per-step attributability | way 1 | the only route with an H artifact and a separate transfer step |
| raw benchmark accuracy of the imputation | way 2 | end-to-end training beats translation tacked onto a fixed embedding |
| scale, with no custom code | way 3 or 5 | a sparse matmul or an aggregation — both are linear, neither needs OT |
| a defensible network rather than a matrix | way 4 | the ATAC evidence stays in the model and regulons become the unit |
| you do not trust imputed values at all | way 6 or 4 | neither one ever produces imputed gene expression |

They compose. Ways 3, 5 and 6 can be stacked: gene activity for the ATAC rows, aggregated into metacells, with TF activity as a second regulator view. That combination needs no optimal transport, no deep model and no custom transfer, which is what makes it the natural candidate at 10⁶ cells.
Way 1 is the only one where every step can be graded independently — which is the reason this deck describes it in detail. That is a design trade, not a claim that it is the most accurate.

**Speaker notes**

> The decision rule is to name the binding constraint first. If it is attributability, way 1 is the only option. If it is raw imputation accuracy, way 2 wins on benchmarks. If it is scale with no custom code, ways 3 and 5 are linear and need no optimal transport. If the deliverable is a defensible network rather than a matrix, way 4. And the routes compose: gene activity for the ATAC rows, aggregated into metacells, with TF activity as a second regulator view, needs no OT, no deep model and no custom transfer. That combination is the natural candidate at a million cells, and it is worth saying out loud that way 1 is not the most accurate route - it is the most inspectable one.

---

## Slide 30 — Validation of the pipeline itself

**On the slide**

scMultiSim lets us grade all three steps against one known answer

**known network** — we specify it

**scMultiSim** — RNA + ATAC + noise

**our four steps** — integration → GRN

**recovered network** — compared

| Step | Metric | Truth |
|---|---|---|
| 1  Integration | FOSCTTM | the simulated pairing |
| 2  Imputation | MSE, pFCS | the expression that was hidden |
| 3  GRN | AUPRC, early precision, hub recovery | the network that was planted |

If the pipeline cannot recover a network it was given — with known structure and controlled noise — then a poor result on real data says nothing about the biology.
scMultiSim: Nat Methods 2025, 10.1038/s41592-025-02651-0.   For plausibility on a finished network: CellOracle in-silico perturbation, Nature 2023.

**Speaker notes**

> Simulation is the only place where we can ask whether the pipeline recovers a network that is definitely there. Answer that before spending compute on real data.

---

## Slide 31 — How the steps fit together

**On the slide**

Each gate prevents one specific failure

| Gate | Protects | From what |
|---|---|---|
| 1 — FOSCTTM + scIB | step 2 | a lookup space where 'nearest neighbour' means nothing, so the average is meaningless |
| 2 — pFCS + MSE + ρ | step 3 | a matrix that is confidently wrong: guessed rows are structured, so the network is structured and wrong |
| 3 — stability + lineage + AUPRC | step 4 | validating noise |
| end-to-end — simulation | the design | building four steps that were never shown to work together |

Four scripts in a row is a sequence.
A pipeline is four steps where each output has been shown good enough to be the next step's input. The gates are the difference — and they are also why the slow step stops being the debugging tool.

**Speaker notes**

> Four scripts in a row is a sequence. A pipeline is four steps where each output has been shown good enough to be the next step's input.

---

## Slide 32 — Run order and cost

**On the slide**

Everything except step 3 runs on CPU in hours

| # | Work | Cost |
|---|---|---|
| 1 | simulated data; grade all three steps | days, CPU |
| 2 | step 1 on real data; run gate 1 | hours |
| 3 | step 2 on real data; run gate 2 | hours to a day |
| 4 | compare two or three imputers, pick one | hours to days |
| 5 | step 3, once: three seeds, per cell type | the expensive part |
| 6 | validation against cell-type-matched truth | days |

**SEQUENCING** — step 3 runs once, after its input has passed a gate — not as the instrument used to test the earlier steps

**Speaker notes**

> Step 3 is the only expensive part, and in this order it runs once, after its input has passed a gate.

---

## Slide 33 — What good looks like

**On the slide**

Set before running anything, so the result cannot be reinterpreted

| Metric | A good result |
|---|---|
| FOSCTTM | 0.14–0.22 (published range); 0.5 is random; 0 is perfect |
| gene-wise ρ, ATAC→RNA | 0.155–0.22 — not 0.9, which suggests leakage |
| pFCS | high, and clearly above the shuffled-reference control |
| MSE | clearly below the same control |
| densification | guessed rows no denser than measured rows |
| network stability | reported differences exceed seed-to-seed noise |
| AUPRC / early precision | above the negative control, on cell-type-matched truth |
| lineage regulators | present in the top edges of every cell-type network |

_[chart: FOSCTTM — lower is better]_

_[chart: ATAC→RNA gene-wise Pearson — published ceiling]_

**Speaker notes**

> Targets are set before running anything, so the result cannot be reinterpreted afterwards.

---

## Slide 34 — Summary

**On the slide**

**1  INTEGRATION**: graded by FOSCTTM and the scIB label/batch metrics; needs paired cells, which multiome data provides.

**2  IMPUTATION**: graded by pFCS above all, then MSE and gene-wise ρ, against a shuffled-reference control.

**3  GRN**: graded by AUPRC and early precision on cell-type-matched truth, with seed stability reported first.

**4  VALIDATION**: the whole pipeline graded once on simulated data with a known network, before real biology.

For every step here there is a published, calibrated way to measure it — three of the four are cheap, because multiome data hands us the answer.

**Speaker notes**

> One line to take away: for every step there is a published, calibrated way to measure it, and three of the four are cheap because multiome data hands us the answer.
