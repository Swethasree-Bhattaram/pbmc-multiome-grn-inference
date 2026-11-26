# MultiVI: the workflow, part by part, and what depends on what

A walk-through of `workflow.png`: what each piece does in plain words, how the pieces depend on each
other, then the equations. Checked against the paper (Ashuach et al., *Nature Methods* 2023) and against
the implementation in **scvi-tools 1.4.2** (`scvi/module/_multivae.py`); where they disagree, §10 says so.

Reading order: §1-§5 for the pieces, **§6 for the dependencies** (including "is the decoder dependent on the
encoder?"), §7 for the math, §8 for the matrix inventory, §9 for one cell end to end.

---

## 0. Your figure, transcribed

```
        Nr x G                                              G x Na
   ┌───────────────┐                                    ┌──────────────────────┐
   │ RNA counts    │                                    │ Imputed matrix on    │
   └──────┬────────┘                                    │ ATAC rows            │
          │                                             └──────────┬───────────┘
   ┌──────▼────────┐         ┌─────────────────────┐              │  (Nr + Na) x G
   │ RNA encoder   │────────▶│ Posterior fusion    │              ▼
   └───────────────┘  Nr x B │ Precision weighted  │      ┌──────────────────────┐
                             │ average             │      │ Imputed matrix on    │
   ┌───────────────┐         └──────────┬──────────┘      │ ATAC rows            │
   │ ATAC encoder  │────────────────────┘                 └──────────┬───────────┘
   └──────▲────────┘  Na x B         ┌───▼──────────┐                │
          │                         │ Latent Space H│                ▼
   ┌──────┴────────┐                └───┬───────┬───┘         ┌──────────────┐
   │ ATAC counts   │                    │       │             │ GRN          │
   └───────────────┘            ┌───────▼──┐ ┌──▼────────┐    └──────────────┘
        Na x P                  │RNA decoder│ │ATAC decoder│
                                └───────┬──┘ └────────────┘
                                        └──────▶ (imputed matrix)
```

### The letters

| figure | meaning | code name | here |
|---|---|---|---|
| `Nr`, `Na` | number of RNA cells, ATAC cells | `n_obs` | `N` |
| `G` | genes | `n_input_genes` | `G` |
| `P` | peaks | `n_input_regions` | `P` |
| `B` | latent width — length of one cell's description | `n_latent` | `L` |
| `H` | the latent matrix, one row per cell | `qz_m` | `z` |

Your figure uses `B` for latent width; I use `L` so it is not confused with minibatch size. Both are 30 in
your runs. `Nr x B` under the RNA encoder is the `mu` it emits.

---

## 1. The task, in plain words

Two measurements of the same tissue:

- **RNA**: per cell, how many reads came off each of 36,601 genes.
- **ATAC**: per cell, how many fragments landed in each of 89,796 open regions.

Some cells were assayed both ways, most only one way. Two problems: the units are incomparable ("UMIs of
GAPDH" vs "fragments in peak 41,203"), and most cells have only one view, so you cannot match rows.

MultiVI's answer: learn a **shared 30-number description per cell** such that either measurement can be
predicted back from it. Cells with RNA get the description from RNA; cells with ATAC get it from ATAC; cells
with both get a blend. Then "impute the missing RNA" is just "run the RNA predictor on that description".

---

## 2. The encoder, in plain words

**In:** one cell — 36,601 counts (RNA) or 89,796 (ATAC), ~95% zeros.
**Out:** 30 numbers `mu` (best guess) and 30 numbers `sigma2` (uncertainty).

Two matrix multiplies and two small heads:

```
36,601  ->  128  ->  128  ->  30 (mu)  +  30 (sigma2)
```

**Where the matrices come from.** Each arrow is a weight table `W` plus a bias. `W1` is 128 x 36,601 = 4.7M
numbers, random at the start, learned by gradient descent. Nothing hand-designed.

**What it learns.** Each row of `W1` is a direction in gene space — a pattern of "these genes move
together". The encoder discovers which patterns describe a cell, combines them in layer two ("this looks
like a B cell"), and outputs 30 coordinates.

**Why `sigma2`.** The 30 numbers are an estimate, not a fact. A shallow cell deserves a vaguer description,
and the model needs a way to say so.

---

## 3. The sample, in plain words

```
z = mu + sigma * noise
```

Jitter the guess by its own uncertainty. Training the decoder on jittered positions stops it memorising
exact cells, which forces the 30 numbers to form a smooth map. At inference the jitter is off: `z = mu`.

---

## 4. The decoder, in plain words

**In:** the 30 numbers. **Out:** 36,601 (RNA) or 89,796 (ATAC) numbers — **not counts**.

The decoder emits **distribution parameters**: for each gene an expected rate, a dispersion, and a dropout
logit; for each peak a probability of being open. The loss then asks "how likely were the observed counts
under that distribution?".

Why: counts are noisy. A gene with true rate 3 reads as 1, 5, or 0. The well-posed target is the rate.

The consequence you will see everywhere: a mean is never sparse. Applied to your 11,000 ATAC cells the
decoder produced a matrix **22.5x denser** than a measured row. Not a threshold — that is what a conditional
mean looks like. It is also why the GRN sits at the noise floor (§9).

---

## 5. Why there are two of each

- **Two encoders, two decoders** because RNA is over-dispersed counts and ATAC is essentially binary; one
  noise model cannot serve both.
- **Masks** let one model hold both. `mask = (rowsum of that modality > 0)`, computed, not supplied. A cell
  with no ATAC counts contributes zero ATAC reconstruction error — it is never asked to explain peaks it
  does not have.
- That asymmetry is the imputation: an ATAC-only cell gets `z` from ATAC alone, and the RNA decoder can still
  be run on that `z`.

---

## 6. How the pieces depend on each other

There are **three different kinds of dependency** here, and mixing them up is the main source of confusion.

1. **Data flow** — at run time, which value is computed from which.
2. **Parameter dependence** — whose weights were *learned* from whose outputs.
3. **Gradient flow** — during training, whose error reaches whose weights.

### 6.1 The graph

Data flows left to right, gradients right to left:

```
  RNA counts ──▶ RNA encoder ──▶ mu_R, sigma_R ──┐
                                                 ├─▶ fusion ──▶ z ──┬──▶ RNA decoder ──▶ px_scale
  ATAC counts ─▶ ATAC encoder ─▶ mu_A, sigma_A ──┘                  │                      │
                                                                    │                      ▼
                                                                    └──▶ ATAC decoder ─▶ p   imputed matrix
                                                                                            (G x Na)
                                          loss ◀── data + masks ◀────────────────────────────┘
                                            │
                                            ├──▶ gradients back through every box above
                                            └──▶ Jeffreys compares mu_R with mu_A directly
```

Two arrows are worth staring at, because they are the non-obvious ones:

- **The Jeffreys term connects the two encoders directly**, not through `z`. It compares `mu_R`/`sigma_R`
  with `mu_A`/`sigma_A` for the same cell. That is the only thing tying the two halves of the model together.
- **The loss touches both decoders and both encoders**, so a change in the ATAC decoder's weights changes
  the RNA encoder's weights by the next gradient step (through `z`). Everything is coupled during training,
  even where it is not coupled at run time.

### 6.2 What each piece needs

| piece | needs at run time | its weights were learned from | receives gradient via |
|---|---|---|---|
| RNA encoder | RNA counts | RNA decoder + KL + Jeffreys | `recon_R` -> `z` -> `mu_R` |
| ATAC encoder | ATAC counts | ATAC decoder + KL + Jeffreys | `recon_A` -> `z` -> `mu_A`, plus Jeffreys |
| fusion | `mu_R`, `sigma_R`, `mu_A`, `sigma_A`, masks, `w` | only `w` is learned (2 or `N x 2` numbers) | `z` into both decoders |
| `z` | the fused `mu`, `sigma2` | nothing — it is computed | carried by both reconstructions |
| RNA decoder | `z` **and** the RNA library size | learned to invert the RNA encoder | `recon_R` + KL |
| ATAC decoder | `z` **and** the ATAC library size | learned to invert the ATAC encoder | `recon_A` |

Note the column that surprises people: **the decoders take two inputs, not one.** `z` says what the cell is;
the library size says how deeply it was sequenced. The RNA rate is `exp(library) * px_scale`, so depth and
biology are separable. That is why the library-size encoders exist even though they produce a single number.

### 6.3 Is the decoder dependent on the encoder?

The honest answer is three answers, because it depends which sense you mean:

| sense | dependent? | why |
|---|---|---|
| **parameters** | **no** | `p(x | z)` is a fixed function with its own weights. It never reads the encoder's weights. You can feed it any 30-vector: a random one, a cell average, an interpolation between two cells, a perturbed one. It will answer. |
| **training** | **yes, circularly** | during training the decoder only ever saw `z` produced by the encoder (plus noise). Its weights are tuned to *that distribution* of `z`. |
| **run time** | **sequential, not dependent** | encoder -> `z` -> decoder is a pipeline, but each step is a separate function. The decoder cannot start until the encoder finishes; it does not need the encoder to exist as a parameter set. |

The practical consequences, in order of how likely they are to bite you:

1. **You can decode anything.** `z = 0` (the prior mean), the midpoint of two cells, a path between cell
   types. The decoder is a standalone generator. This is what makes perturbation and interpolation analyses
   possible.
2. **You cannot mix a decoder with a foreign encoder.** This is the one that matters for your repo. Your
   pipeline saves `joint_embedding_H.npy` for both methods, both `N x L`, both called `H`. They are **not the
   same space**. A MultiVI decoder trained on MultiVI's `z` has never seen a scSAGA `H` and will produce
   nonsense on it — the numbers have the same shape and no shared meaning. Never cross them, and never
   compare them directly beyond scale-free metrics.
3. **The decoder degrades gracefully with distance.** Because it was trained on the encoder's `z`
   distribution, decoding a `z` from a *different* dataset (reference mapping) works only if that dataset's
   cells land in the same region. That is exactly what the Jeffreys term buys you, and why the unpaired case
   fails: without it, the ATAC encoder's `z` for the same cell sits far from the RNA encoder's, and the RNA
   decoder produces a blur (measured FOSCTTM 0.43 vs 0.003).

In one sentence: **the decoder does not depend on the encoder, but it depends on the *space* the encoder
produces — and that space is a property of the pair, not of either alone.**

### 6.4 What breaks when a piece is missing

| missing piece | consequence | measured here |
|---|---|---|
| co-assayed cells | Jeffreys = 0, the two encoders drift, one latent space becomes two | FOSCTTM 0.43 (random 0.5) vs 0.003 paired |
| RNA for a cell | `z` from ATAC only; the RNA decoder's output is a prediction, not a reconstruction | the imputation path |
| ATAC for a cell | `z` from RNA only; the ATAC mask silences its peak error | — |
| KL term | no pressure toward a common space; every cell gets a private coordinate | not run |
| correct library size | rate scaled by the wrong depth; gene means shift uniformly | gene mean shift 0.336 sd (MultiVI) |
| matching latent space | decoder emits a blur or noise | see 6.3, point 2 |

---

## 7. The mathematics

Notation: modality `m` in `{R, A}`; cell `i`; `x_i^R` gene counts (length `G`); `x_i^A` peak counts
(length `P`); `L = 30`; `H = 128`; minibatch size `B`.

### 7.1 Encoder (per modality, separate weights)

```
h1 = LayerNorm( LeakyReLU( W1 x + b1 ) )      W1: H x G (or H x P)        B x H
h2 = LayerNorm( LeakyReLU( W2 h1 + b2 ) )     W2: H x H                  B x H
mu     = W_mu h2 + b_mu                       W_mu: L x H                B x L
sigma2 = exp( W_v h2 + b_v )                  W_v:  L x H                B x L
```

From `scvi.nn.Encoder` (used at `_multivae.py:367`): activation is **LeakyReLU**; **LayerNorm** after each
layer; dropout 0.1; variance is `exp(linear)` so it is positive by construction. The code stores a
**variance**, not a standard deviation.

A second encoder per modality produces one number, the library size:

```
l = softplus( w_l h' + b )                    length 1 per cell
```

### 7.2 Sampling

```
eps ~ Normal(0, I)              B x L
z   = mu + sqrt(sigma2) * eps   B x L
```

Randomness sits outside the network, so gradients reach `mu` and `sigma2` normally.

### 7.3 Fusion

```
W_m = masked_softmax(w_m, mask_m)      across modalities, sum to 1, absent modality -> 0
mu  = sum_m W_m * mu_m                 B x L
var = sum_m sqrt(W_m) * var_m          B x L
```

`w` is per-modality: `equal` (0.5/0.5, default), `universal` (2 learned numbers), or `cell` (an `N x 2`
learned table).

A cell with one modality keeps that belief exactly (`W = 1` for the present one). This is **not** precision
weighting and **not** a product of experts — nothing divides by a variance. The paper writes it the same way:
`Z_c = w_A Z_c^A + w_R Z_c^R`, `sum_m w_m = 1`. The variance line is a heuristic pool (`sqrt` weights), exact
only for a single modality.

### 7.4 RNA decoder

```
h        = LeakyReLU stack(z)                  B x H    (two layers, LayerNorm)
px_scale = softmax( W_s h + b_s )              B x G    proportions over genes, sum to 1
rate     = exp(l) * px_scale                   B x G    expected counts
theta    = exp( W_t h + b_t )                  G        per-gene inverse dispersion
pi       = W_d h + b_d                         B x G    dropout logits
```

Loss per cell, summed over genes:

```
ZINB (default):  -log p(x | rate, theta, pi) = -sum_g log ZINB(x_g ; rate_g, theta_g, pi_g)
NB:              drop pi
Poisson:         rate only
```

The `softmax` is why `px_scale` is a **within-cell proportion**: it says what share of that cell's reads went
to each gene; multiplying by `exp(l)` makes expected counts. `log1p(1e4 * px_scale)` is the imputed
expression your pipeline stores.

### 7.5 ATAC decoder

```
h = LeakyReLU stack(z)                 B x H
p = sigmoid( W_p h + b_p )             B x P     probability each peak is open
r = sigmoid(region_factors)            P         learned per-peak accessibility prior
```

Loss per cell, summed over peaks, against binarised counts:

```
BCE = -sum_j [ x_j log(p_j l r_j) + (1 - x_j) log(1 - p_j l r_j) ] ,   x_j = 1[count > 0]
```

### 7.6 Objective

Per minibatch, averaged over cells:

```
L = mean_i [  mask_R(i) * (-log p(x_i^R | z_i))            RNA reconstruction
            + mask_A(i) * BCE(x_i^A, p_i)                  ATAC reconstruction
            + beta * KL( q(z_i|x_i) || Normal(0,I) )       keep the map tidy
            + Jeffreys( q_R(z_i) || q_A(z_i) ) ]           make the encoders agree
```

`beta = 1`; no separate weight on Jeffreys in scvi-tools 1.4.2 — added to the KL directly (`_multivae.py:868`).

```
KL = 0.5 * sum_k [ mu_k^2 + sigma_k^2 - log(sigma_k^2) - 1 ]

J  = KL(q_R||q_A) + KL(q_A||q_R)
   = 0.5 * sum_k [ sigma_R^2/sigma_A^2 + sigma_A^2/sigma_R^2 - 2
                   + (mu_R - mu_A)^2 (1/sigma_R^2 + 1/sigma_A^2) ]
```

`J` is computed only where `mask_R & mask_A` (`_multivae.py:934-970`). With no co-assayed cells it is exactly
zero — the failure mode in §6.4.

### 7.7 Optimisation

Shuffled minibatches of 128, Adam, 200 epochs, early stopping on held-out ELBO. One pass over 11,000 cells is
86 gradient steps. Nothing is closed-form; every matrix in §8 is updated by gradient descent.

---

## 8. Where every matrix is, and how big

`G = 36,601`, `P = 89,796`, `H = 128`, `L = 30`. Counted by instantiating the module.

| module | layers | parameters |
|---|---|---|
| `z_encoder_expression` | 36,601 -> 128, 128 -> 128, 128 -> 30 (mean), 128 -> 30 (var) | 4,709,308 |
| `l_encoder_expression` | 36,601 -> 128, 128 -> 128, 128 -> 1 | 4,701,697 |
| `z_decoder_expression` | 30 -> 128, 128 -> 128, three heads 128 -> 36,601 | 14,185,067 |
| `z_encoder_accessibility` | 89,796 -> 128, 128 -> 128, 128 -> 30 x 2 | 11,518,268 |
| `l_encoder_accessibility` | 89,796 -> 128, 128 -> 128, 128 -> 1 | 11,510,657 |
| `z_decoder_accessibility` | 30 -> 128, 128 -> 128, 128 -> 89,796 | 11,604,164 |
| protein stubs (0 proteins) | — | 85,948 |
| **total** | | **58,441,506** |

The ATAC encoder's first layer alone is 11.5M numbers — 20% of the model in one matrix, because it must map
89,796 peaks to 128 numbers. RNA pays for three decoder heads, ATAC for one, which is why 36,601 genes cost
14.2M and 89,796 peaks cost 11.6M.

---

## 9. One cell, end to end

An ATAC-only cell from your 11,000:

1. **Encode.** 89,796 peak counts (~2,500 non-zero) -> `W1` (128 x 89,796, ~11.5M multiply-adds) -> two more
   small layers -> `mu` (30) + `sigma2` (30). Its RNA mask is 0, so the RNA encoder contributes nothing.
2. **Fuse.** One belief: `W_A = 1`, `W_R = 0`. `z = mu`.
3. **Decode.** `z` -> RNA decoder (30 -> 128 -> 128 -> 36,601, three heads, ~14M multiply-adds) -> `px_scale`,
   a 36,601-long vector of proportions summing to 1.
4. **Export.** `log1p(1e4 * px_scale)` becomes one column of the `36,601 x 11,000` imputed matrix, stacked
   under the measured RNA rows to make `(Nr + Na) x G` for Arboreto.

Why the GRN result follows: every entry is a conditional mean from a softmax — no zeros, and neighbouring
cells (nearby `z`) get nearly identical rows. GRNBoost2 looks for co-variation across cells, and a smooth
function of a 30-dimensional latent has strong smooth co-variation that is not in the data. That is the
mechanism behind AUPRC 0.006 against prevalence 0.005 and early precision 0.000 — not a broken network code,
but an input that is too regular.

---

## 10. Three places the figure and the code disagree

1. **"Posterior fusion: Precision weighted average"** -> the code does a **masked weighted mean**
   (`masked_softmax`, then a weighted sum; variances weighted by `sqrt(w)`). No `1/sigma^2` anywhere. The paper
   agrees with the code. Relabel the box.
2. **The arrow from the ATAC decoder into the imputed expression matrix** -> the `G x Na` matrix comes from the
   **RNA decoder** applied to ATAC cells' `z`. The ATAC decoder produces a `P x Na` peak matrix your pipeline
   never uses. Remove that arrow.
3. **`n_hidden` default** is `min(128, sqrt(n_features))`, so it is 128 for any realistic dataset; it does not
   grow with the peak count.

One more, not in the figure but affecting results: the **default likelihood is ZINB**, so the dropout head is
used. Under `gene_likelihood="nb"` that head is allocated and ignored — a 4.7M-parameter no-op.

---

## 11. Glossary

| term | plain meaning |
|---|---|
| encoder | funnel from a cell's counts to 30 numbers + 30 variances |
| decoder | reverse funnel from 30 numbers to distribution parameters for every feature |
| latent `z` | the 30 numbers; the cell's coordinates |
| `mu`, `sigma2` | per-coordinate guess and variance from one encoder |
| mask | `rowsum(counts) > 0`; which modality a cell actually has |
| fusion | combining the two beliefs into one `z` |
| library size `l` | sequencing depth of a cell; the decoder's second input |
| `px_scale` | per-cell proportions over genes; the imputed expression |
| `rate` | `px_scale` x library size; expected counts |
| dispersion `theta` | over-dispersion of a gene's counts |
| dropout `pi` | ZINB zero-inflation logit |
| `p` | per-peak probability of being open |
| region factor | learned per-peak accessibility prior |
| KL | penalty keeping the map centred and smooth |
| Jeffreys | penalty forcing the two encoders to agree on co-assayed cells |
| ELBO | the objective: data likelihood minus the penalties |
