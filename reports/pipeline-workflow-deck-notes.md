# Deck notes: from multiome data to a gene regulatory network

Companion to `reports/pipeline-workflow-design.pptx` (20 slides) and to
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

## Slide 9 — Step 2 — imputation methods that use H

**On the slide**

Only this family keeps step 1 and step 2 separate — and it is what you already do

| Method | What it consumes | How it makes gene values | Citation |
|---|---|---|---|
| reverse-imputeKNN  (this repo) | H + the measured RNA + a reference block | 20 nearest reference cells, weights exp(−distance) | run.py |
| LIGER  imputeKNN | H (its H.norm) + the measured RNA | 20 nearest reference cells, weighted average | Welch et al., Cell 2019 |
| moscot | H passed as the shared embedding (joint_attr) | optimal-transport plan between cells; expression transported along it | Klein et al., Nature 2025 |
| Symphony  mapQuery | H + cluster centroids of the reference | maps query cells into H, then nearest-cell transfer | Kang et al., Nat Commun 2021 |
| MOFA+  project_data | H as factors + their gene loadings | projects cells into the factor space, then reconstructs | Argelaguet et al., Genome Biol 2020 |

H + MEASURED RNA  H plus the measured RNA is all these need — which is why the integration step earns its place: without H there is no shared space in which to look up a neighbour

None of these can work from H alone: coordinates do not contain expression. Every one of them is H plus one companion object.

**Speaker notes**

> This is the family your design needs, and it is the one the repo is in. Every entry consumes H plus one companion object, because coordinates alone cannot produce expression values. Our reverse-imputeKNN takes H, the measured RNA and a frozen reference set, finds the twenty nearest reference cells and averages with exp minus distance weights. LIGER's imputeKNN is the same algorithm - it refuses to run if H.norm is missing, then takes the twenty nearest reference cells in that space. moscot accepts an embedding as its joint_attr and solves an optimal-transport plan between cells, then moves expression along that plan; that is the OT analogue of the same idea. Symphony's mapQuery projects query cells into a reference embedding and you transfer from the nearest reference cells. MOFA+ projects into an existing factor space and reconstructs. The point for the talk: the integration step is load bearing here. Without H there is no shared space in which a neighbour lookup means anything, so step 1 is not decoration - it is the only reason the imputation works.

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

## Slide 16 — Validation of the pipeline itself

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

## Slide 17 — How the steps fit together

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

## Slide 18 — Run order and cost

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

## Slide 19 — What good looks like

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

## Slide 20 — Summary

**On the slide**

**1  INTEGRATION**: graded by FOSCTTM and the scIB label/batch metrics; needs paired cells, which multiome data provides.

**2  IMPUTATION**: graded by pFCS above all, then MSE and gene-wise ρ, against a shuffled-reference control.

**3  GRN**: graded by AUPRC and early precision on cell-type-matched truth, with seed stability reported first.

**4  VALIDATION**: the whole pipeline graded once on simulated data with a known network, before real biology.

For every step here there is a published, calibrated way to measure it — three of the four are cheap, because multiome data hands us the answer.

**Speaker notes**

> One line to take away: for every step there is a published, calibrated way to measure it, and three of the four are cheap because multiome data hands us the answer.
