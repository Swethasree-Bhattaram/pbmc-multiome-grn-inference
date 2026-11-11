# MultiVI: architecture, time and space complexity, and scaling to 10^6 cells

Scope: `scvi.module.MULTIVAE` as shipped in **scvi-tools 1.4.2** (torch 2.15.0), the model behind
`scripts/stepwise_eval.py`. Every parameter count below was produced by instantiating the module and
counting `p.numel()`; the timing anchor is the benchmark in §8. Anything extrapolated rather than
measured is labelled as such.

## 1. Notation

| symbol | meaning | repo value |
|---|---|---|
| `N` | cells | 8,381 RNA / 8,728 ATAC (17,109 collapsed rows) |
| `G` | genes (RNA features) | 36,601 |
| `P` | ATAC regions / peaks | 89,796 |
| `h` | hidden width per layer (`n_hidden`) | 128 |
| `L` | latent width (`n_latent`) | 30 |
| `B` | minibatch size | 128 (scvi default) |
| `X_r` | RNA counts, sparse `N x G` | |
| `X_a` | ATAC counts, sparse `N x P`, binarised in the loss | |
| `m_r`, `m_a` | modality masks, `N` each, `m = (rowsum X > 0)` | |
| `b` | batch label per cell | none (single experiment) |

The masks come from the data, not from the caller: a cell with zero counts in a modality is treated as
that modality missing. This is the mechanism that lets one model hold RNA-only, ATAC-only and paired
cells.

## 2. The model in six steps

### Step 0 — inputs

```
X_r   N x G    sparse raw RNA counts
X_a   N x P    sparse ATAC counts
b     N        batch label  -> one-hot covariate concatenated to every encoder input
m_r   N        RNA present
m_a   N        ATAC present
```

### Step 1 — one encoder per modality, separate weights

```
x_r  B x G -> FC(G->h) -> ReLU -> FC(h->h) -> ReLU      B x h
           -> Linear(h->L), Linear(h->L)               mu_r, logvar_r    B x L
x_r  B x G -> library-size encoder (G->h->1, softplus) l_r               B x 1

x_a  B x P -> FC(P->h) -> ReLU -> FC(h->h) -> ReLU      B x h
           -> Linear(h->L), Linear(h->L)               mu_a, logvar_a    B x L
x_a  B x P -> library-size encoder (P->h->1)            l_a               B x 1
```

Each encoder sees **only its own modality**. The concatenated `G+P` vector never exists as an
activation; the two modalities meet for the first time in step 2. Consequence: an RNA-only cell costs
nothing on the ATAC encoder, and vice versa.

Measured cost of this split (repo widths): RNA encoders 9.41 M params, ATAC encoders 23.03 M params.

### Step 2 — fuse the two posteriors into one joint `q(z)`

This is a **weighted average**, not a product of experts, and not precision weighting. Verified in
`scvi/module/_multivae.py` (`mix_modalities`, called from `inference`) and consistent with the paper,
which writes `Z_c = w_A Z_c^A + w_R Z_c^R` with `sum_m w_m = 1`.

```
mu_r, var_r, mu_a, var_a, m_r, m_a, w
      |
      +--> weights = masked_softmax(w, m)     # zero the absent modality, renormalise to sum 1
      |
mu   = sum_m  weights_m * mu_m                                            B x L
var  = sum_m  sqrt(weights_m) * var_m                                     B x L
```

The variance line really is `sqrt(weights)` then `sum`, not `sum(weights * var)` and not
`sqrt(sum(weights * var))`; it is a heuristic pool rather than a Gaussian product, and it is exact only
when a cell has one modality (there the single weight is 1 and `var = var_m`). Two consequences worth
knowing:

- **`w` is a mixture weight, not a precision.** Nothing here divides by a variance, so there is no
  `sigma -> 0` blow-up and no per-cell normaliser to stabilise. The fusion is cheap and stable.
- **Cost is elementwise on `B x L`** — about ten kernels over 128 x 30 floats. Measured against the
  encoder's `B x P x h` GEMM (128 x 250,000 x 128), the fusion is ~5 orders of magnitude cheaper. It is
  the conceptual core of the model and the place the two modalities meet, but it is not a compute cost.

`w` is one of:

| `modality_weights` | shape | cost |
|---|---|---|
| `equal` (default) | scalar | 0 |
| `universal` | `2` | 2 |
| `cell` | `N x 2` | 8 MB at 1M cells (fp32) |

A cell with one modality keeps that encoder's belief, because `masked_softmax` gives the absent modality
weight 0 and renormalises the present one to 1.

The **training-time** counterpart of this step is the Jeffreys penalty in the loss, and that one is *not*
cheap: it forms a `B x B x L` object per modality pair (see §6). That is the term that caps batch size.

### Step 3 — reparameterised sample

```
z = mu + sqrt(var) * eps ,   eps ~ N(0, I)        B x L
```

### Step 4 — decode both modalities from the same `z`

```
z  B x L -> FC(L->h) -> ReLU -> FC(h->h) -> ReLU        B x h
   RNA head:  FC(h->G) -> normalised mean (px_scale)     B x G
              FC(h->G) -> dispersion (px_r)              G
              FC(h->G) -> dropout logits (px_dropout)    B x G
   ATAC head: FC(h->P) -> Bernoulli p                    B x P
              * region_factors (P)  * library size l_a
```

The RNA decoder allocates **three** `h x G` matrices (mean, dispersion, dropout) even when the
likelihood is plain NB, so the dropout head is paid for but unused. That is why the RNA decoder is
14.19 M params against the ATAC decoder's 11.60 M at the same width.

### Step 5 — loss

```
L = - sum_genes  m_r * log NB(x_r  | px_scale, px_r, px_dropout, l_r)
    - sum_peaks  m_a * BCE(x_a>0 | p)
    + beta * KL( q(z|x) || N(0,I) )
    + lambda * sum_{cells with m_r AND m_a}  Jeffreys( q_r(z) || q_a(z) )
```

The last term is what makes MultiVI multi-omic rather than two independent VAEs. It acts only on
co-assayed cells: with no cell carrying both modalities it is identically zero, and the two encoders
drift apart. Measured on synthetic paired data: FOSCTTM 0.003 with collapsed paired rows vs 0.43 when
RNA and ATAC sit on separate rows (`reports/stepwise-pbmc10k-subset.md`).

Its cost is `O(B^2)` per batch (a batch x batch divergence), not `O(N^2)`: 16,384 entries at `B = 128`.

### Step 6 — minibatch Adam

`max_epochs=200`, early stopping on validation ELBO, `batch_size` left at the scvi default.

```
loss.backward() -> Adam step        per batch, O(N/B) batches per epoch
```

### Export

```
get_latent_representation()  ->  mu_z           N x L
decode ATAC rows             ->  px_scale       n_atac x G   (streamed, batch by batch)
```

## 3. Parameter count

Closed form, ignoring `O(h^2 + hL)` terms:

```
params  =  5*h*G  +  3*h*P
           |          |
           |          +-- ATAC: z-encoder, library encoder, Bernoulli decoder
           +-------------- RNA: z-encoder, library encoder, 3 decoder heads
```

Measured, module by module (scvi-tools 1.4.2, `n_layers_encoder = n_layers_decoder = 2`):

| module | repo config (P=89,796) | atlas (P=250,000, h=128) | atlas (P=250,000, h=500 default) |
|---|---|---|---|
| `z_encoder_expression` | 4,709,308 | 4,709,308 | 18,581,560 |
| `l_encoder_expression` | 4,701,697 | 4,701,697 | 18,552,001 |
| `z_decoder_expression` | 14,185,067 | 14,185,067 | 55,277,303 |
| `z_encoder_accessibility` | 11,518,268 | 32,024,380 | 125,281,060 |
| `z_decoder_accessibility` | 11,604,164 | 32,270,480 | 125,516,000 |
| `l_encoder_accessibility` | 11,510,657 | 32,016,769 | 125,251,501 |
| protein stubs (unused) | 85,948 | 85,948 | 1,079,560 |
| **total** | **58,441,506** | **120,280,250** | **469,825,586** |
| formula `5hG+3hP` | 57,906,304 | 119,424,640 | 466,502,500 |
| Adam footprint (param + grad + 2 moments, fp32) | 0.94 GB | 1.92 GB | 7.52 GB |

Two things this makes obvious:

1. **The peak axis is the cost driver.** `P` enters three times, `G` five times, and `P` is the larger
   number in any atlas-scale dataset. Halving peaks halves two thirds of the model.
2. **`n_hidden` defaults to `sqrt(n_regions)`.** Leave it unset at 250k peaks and the model is 3.9x
   bigger than `n_hidden=128` (469.8 M vs 120.3 M params). Set `n_hidden` explicitly. The repo already
   does (`--n-hidden 128`).

None of these are large numbers. 58 M params is a small model by 2026 standards; the memory story at
1M cells is the data, not the weights.

## 4. Time complexity

```
FLOPs per cell per epoch  ~  6 * params   (2x forward, 4x forward+backward)
time per epoch            ~  O(N * h * (G + P))        linear in cells, genes and peaks
```

Linear, with no `O(N^2)` or `O(N*N)` term anywhere in the model. Same asymptotics as scVI; unlike
scSAGA / SCOT / Pamona, whose coupling step is `O(N^2)` in cells.

Measured anchor: **1.585 ms per cell per epoch** at the repo's width (`G+P = 126,397`, `h = 128`),
8 CPU threads, 2,000 rows, 3.17 s/epoch. Extrapolated linearly:

| cells | s/epoch (8 CPU threads) | 100 epochs | 200 epochs |
|---|---|---|---|
| 2,000 | 3.17 (measured) | 5 min | 11 min |
| 8,728 (`atac8k`) | 13.8 | 23 min | 46 min |
| 17,109 (collapsed rows) | 27.1 | 45 min | 1.5 h |
| 100,000 | 159 | 4.4 h | 8.8 h |
| 1,000,000 | 1,585 | 44 h | 88 h |

Per fit. The pipeline fits twice (all cells; then with held-out cells excluded), so double these.

GPU, **extrapolated, not measured here**: dense GEMMs of this shape run 30-60x faster on one modern
GPU than on 8 CPU threads, which puts 1M cells at roughly 26-53 s/epoch, i.e. 0.7-1.5 h per 100-epoch
fit. That is an inference from the CPU FLOP rate, not a benchmark. What is measured is the *ratio*
between thread counts, not the GPU number.

## 5. Space complexity

| item | formula | repo config | 1M cells, 250k peaks |
|---|---|---|---|
| model + optimizer | `16 * params` bytes | 0.94 GB | 1.92 GB |
| minibatch activations | `B*(G+P)*4` bytes | 65 MB | 147 MB |
| RNA matrix on disk | `nnz_r * 8` (fp32 data + int32 idx) | — | ~15 GB at 5% nnz |
| ATAC matrix on disk | `nnz_a * 8` | — | ~20 GB at 1% nnz |
| decoded export | `N*G*4` bytes | 1.4 GB at 17k cells | **146 GB** |

Sparsity matters: `nnz = N * F * s`, so bytes `= N * F * s * 8`. At 1M cells and 250k peaks, 1% nnz
is already 2.5e9 non-zeros.

## 6. Bottleneck

Ordered by where the wall actually is, at each scale.

**10^4 - 10^5 cells: compute.** The `h*(G+P)` layers dominate. `P` is the lever; nothing else matters.

**10^6+ cells: the data path, in this order.**

1. **Sparse-to-dense collation.** Every minibatch materialises `B x (G+P)` floats (65 MB at repo width,
   147 MB at 250k peaks) built from a sparse store, on the CPU side. One epoch at 1M cells streams
   ~100 GB of dense rows. This is what Lightning's `does not have many workers` warning is about, and
   it decides whether the accelerator runs or starves. Fix: backed/memmapped AnnData, several
   dataloader workers, larger `batch_size`.
2. **Export.** `N x G` dense is 146 GB at 1M cells and cannot be materialised once. Must be decoded and
   written in row chunks. `scripts/stepwise_eval.py` currently holds `genes x cells` in RAM, which is
   fine at 17k cells (1.4 GB) and not at 1M.
3. **Peak axis.** Linear in `P` for params, FLOPs and I/O. Reducing peaks is a straight 1:1 saving, and
   it must be decided before the fit because `P` is the decoder's output axis.

Not bottlenecks: no coupling matrix, no kNN graph, no per-cell quadratic object. The alignment penalty
is `B^2` inside a batch (16k entries at `B=128`); per-cell modality weights are `2N` floats. Raising
`batch_size` into the thousands is the only way that term becomes visible.

## 7. Does it scale to more than a million cells?

**Architecturally yes.** Nothing in the model is `O(N^2)`; training is minibatch and the export is
streamable. The scVI family has been trained at 10^7+ cells in atlas work.

Two honest caveats specific to MultiVI:

1. **No published MultiVI result above 10^5 cells.** The paper (Ashuach et al., *Nat. Methods* 2023) tops
   out at **89,655 cells** — mixed-source PBMC: 47,148 ATAC-only (Satpathy), 30,495 RNA-only (Ding),
   12,012 paired (10x). The remaining datasets are 10x PBMC 10k, DOGMA-seq and TEA-seq. The word
   "million" does not appear anywhere in the article, and it explicitly lists "datasets with a small
   number of cells" as a limitation (the opposite wall). The ">1 million cells" claim comes from the
   **scvi-tools model page** (`docs.scvi-tools.org`), which is documentation, not a published benchmark —
   the same page also says MultiVI "effectively requires a GPU for fast inference", which matches the
   CPU table in §4. The 10^6 figure in this document is an extrapolation from the architecture.
2. **Co-assayed cells are a hard requirement.** The alignment term only acts where a cell has both
   modalities. A 1M-cell atlas that is 99% single-modality gives the two encoders almost no signal to
   agree on, and the model degenerates toward two independent VAEs — the exact failure measured above
   (FOSCTTM 0.43).

Also note what does *not* scale with it: the downstream steps. FOSCTTM and label-transfer metrics need a
kNN over 1M cells, and the GRN step is one regression per gene column. Those are the expensive parts of
this pipeline at 10^6, not the VAE.

## 8. Reproducing the numbers

Parameter counts (§2, §3):

```python
from scvi.module import MULTIVAE
m = MULTIVAE(n_input_genes=36601, n_input_regions=89796,
             n_hidden=128, n_latent=30,
             n_layers_encoder=2, n_layers_decoder=2)
print(sum(p.numel() for p in m.parameters()))          # 58,441,506
for n, mod in m.named_children():
    print(n, sum(p.numel() for p in mod.parameters()))
```

Timing anchor (§4) — 2,000 rows at repo width, threaded vs single-threaded:

```
OMP_NUM_THREADS=8   -> 3.17 s/epoch
OMP_NUM_THREADS=1   -> 7.14 s/epoch        (2.3x penalty, hence MULTIVI_THREADS in run_stepwise.slurm)
```

Run on a Mac (M-series, CPU only). Absolute values do not transfer to PACE Phoenix; the per-cell cost
and the thread ratio do.

## 9. Where this lives in the repo

| file | role |
|---|---|
| `scripts/stepwise_eval.py` | fits MultiVI (`_multivi_fit`), decodes it (`_multivi_decode`) |
| `README_STEPWISE.md` | the experiment and the PACE runbook |
| `reports/stepwise-pbmc10k-subset.md` | measured 450-cell two-method results |
| `docs/multivi-architecture-and-scaling.md` | this document |
