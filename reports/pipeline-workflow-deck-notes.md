# Deck notes: from multiome data to a gene regulatory network

Companion to `reports/pipeline-workflow-design.pptx` (18 slides) and to
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
- paired cells held out

**OUT**
- H: cells × 30 numbers
- row-range map per dataset
- no genes, no imputation

**2  IMPUTATION**

**IN**
- H  (from step 1)
- real RNA, one gene axis
- list of ATAC rows
- held-out cells frozen

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

> Read it as a factory line: each output is the next step's input. Only two objects matter, H and X, and both can be inspected without running a network.

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

**NOTE** — fix the ATAC feature space (peaks or gene activity) before comparing integrators — ATAC integration is more sensitive to it than RNA is

**Speaker notes**

> Integration never sees genes or peaks, only 50 summaries per cell. The paired cells are held out from the start, because they are how this step gets graded.

---

## Slide 5 — Step 1 — Integration: methods

**On the slide**

Six options; choose with the benchmark, not by intuition

| Method | What it does | Choose it when | Citation |
|---|---|---|---|
| scSAGA | optimal-transport alignment to an anchor dataset | current baseline | Aluru lab (scSAGA) |
| scGT | graph transformer; paired + unpaired; label transfer | scaling to millions of cells | Bioinformatics 2025 |
| scMODAL | aligns using a few known linked features | few feature links available | Nat Commun 2025 |
| SuperMap | learns the cross-modal map from unpaired data | ATAC-only datasets are added | PNAS 2026 |
| scHPGT | regulatory-prior attention; resists over-alignment | cell states must stay distinct | Bioinformatics 2026 |
| MultiVI | one probabilistic model for RNA + ATAC | steps 1 and 2 as a single model | Nat Methods 2023 |

**READ FIRST** — Liu et al., Nat Methods 2025 — 40 methods, with a decision guideline for dimension reduction, batch correction, clustering, label transfer and imputation

**Speaker notes**

> Do not choose by intuition: the 2025 Nature Methods benchmark covers these exact tasks and publishes a decision rule. MultiVI is the option if steps 1 and 2 should be one model.

---

## Slide 6 — Step 1 — Integration: how good is it?

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

## Slide 7 — Step 2 — Imputation: why this step is risky

*Averaging is smoothing, and smoothing manufactures correlation*

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

## Slide 8 — Step 2 — Imputation: input and output

**On the slide**

One gene axis, one normalization, one frozen held-out set

|  | Item | Shape / form | Why it matters |
|---|---|---|---|
| IN | H.npy + blocks.json | cells × 30 | the lookup space from step 1 |
| IN | real RNA | cells × genes, ONE gene axis | the values being copied |
| IN | normalization | declared, e.g. log1p CPM/1e4 | measured and guessed rows must share a scale |
| IN | ATAC row list | indices | which rows are being filled |
| IN | held-out paired cells | frozen index | never in the reference, or the test is void |
| OUT | X.npy | cells × genes, float32 | the only input to step 3 |
| OUT | mask.npy | cells × genes, bool | per cell per gene: measured or guessed |
| OUT | genes.txt, norm.txt | one axis, one normalization | written once, never re-derived |
| OUT | frozen reference set | cell indices | keeps the held-out test honest |

**WHY THE MASK** — without the mask nobody downstream can point at a row and say 'that one is an average' — so nobody can ask whether step 2 helped

**Speaker notes**

> Three constraints: one gene axis, one normalization, one frozen held-out set. The mask is the addition that makes this step testable - with it, step 3 can be run with and without the guessed rows.

---

## Slide 9 — Step 2 — Imputation: methods

**On the slide**

The published ceiling for this task is low — expect 0.2, not 0.9

| Method | In plain words | Choose it when | Citation |
|---|---|---|---|
| kNN average | average the 20 nearest RNA cells | simplest; see the cost on the previous slide | this repo |
| scGALA | step-1 alignment plus a learned mapper | closest drop-in for alignment-then-map | Nat Commun 2025 |
| scDCT | conditional diffusion, both directions | topology preservation matters most | Brief Bioinform 2025 |
| Cisformer | pretrained; 10,000 open regions per cell | more than ~30k cells | Genome Biol 2025 |
| scButterfly | two aligned VAEs, paired or unpaired | simplest learned model; good control arm | Nat Commun 2024 |
| MultiVI | one model does steps 1 and 2 | paired data only | Nat Methods 2023 |
| Polarbear | original semi-supervised translator | the number to beat | J Comput Biol 2022 |

**CEILING** — ATAC→RNA gene-wise Pearson ≈ 0.22 (Polarbear), 0.155 (BABEL), on SNARE-seq.   0.9 suggests leakage; 0.05 means the step produced nothing.

**Speaker notes**

> The ceiling is low. Polarbear reports about 0.22 gene-wise Pearson for ATAC to RNA, BABEL about 0.155. If a method reports 0.9, suspect leakage rather than excellence. If it reports 0.05, the step produced nothing. Both are knowable in an afternoon.

---

## Slide 10 — Step 2 — Imputation: how good is it?

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

MSE alone is not enough: a blurred matrix can have low error and still destroy correlation structure.   Definitions follow the imputation task in Liu et al., Nat Methods 2025; scIMC (NAR 2022) and the scTranslation benchmark (ACM SIGKDD 2026) are ready-made harnesses.

**Speaker notes**

> pFCS asks the exact question step 3 depends on: after imputation, do genes still correlate the way they really correlate? A blurred matrix can have low MSE and still destroy correlation structure, and correlation structure is the only thing the network model reads.

---

## Slide 11 — Step 3 — GRN inference: input and output

**On the slide**

The model does not know which rows were guessed

|  | Item | Shape / form | Why it matters |
|---|---|---|---|
| IN | X.npy + mask.npy | cells × genes | from step 2 |
| IN | regulator list | curated TFs; count recorded | what may appear on the left of an edge |
| IN | target rule | stated rule; count recorded | what may appear on the right |
| IN | cell grouping | one network per cell type, or pseudobulk | a mixed network is dominated by cell-type differences |
| OUT | edges.csv | TF, target, importance | one file per cell type |
| OUT | run_info.txt | cells, columns, seeds, software versions | required for any comparison |
| OUT | three seeds | the same run repeated | makes seed-to-seed noise measurable |
| OUT | top-50k cap | ranked edges used for scoring | the standard evaluation range |

**RECORD** — write both counts down — two runs with different target counts produce edge sets that cannot be compared, and the difference is invisible without them

**Speaker notes**

> The booster does not know which rows were guessed. Two more things: write down the regulator and target counts (otherwise two runs cannot be compared), and compute one network per cell type rather than one over a mixture.

---

## Slide 12 — Step 3 — GRN inference: methods

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

## Slide 13 — Step 3 — GRN inference: how good is it?

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

> Use AUPRC rather than AUROC: with a few thousand true edges among hundreds of thousands of possible pairs, AUROC looks flattering while AUPRC tells the truth. Score only the top-ranked edges, and report stability before anything else - a difference smaller than the seed noise is not a measurement.

---

## Slide 14 — Validation of the pipeline itself

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

## Slide 15 — How the steps fit together

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

## Slide 16 — Run order and cost

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

## Slide 17 — What good looks like

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

## Slide 18 — Summary

**On the slide**

**1  INTEGRATION**: graded by FOSCTTM and the scIB label/batch metrics; needs paired cells, which multiome data provides.

**2  IMPUTATION**: graded by pFCS above all, then MSE and gene-wise ρ, against a shuffled-reference control.

**3  GRN**: graded by AUPRC and early precision on cell-type-matched truth, with seed stability reported first.

**4  VALIDATION**: the whole pipeline graded once on simulated data with a known network, before real biology.

For every step here there is a published, calibrated way to measure it — three of the four are cheap, because multiome data hands us the answer.

**Speaker notes**

> One line to take away: for every step there is a published, calibrated way to measure it, and three of the four are cheap because multiome data hands us the answer.
