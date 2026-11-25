# MultiVI: what the encoder and decoder actually do

A walk-through of your `workflow.png`, in plain words first, then the equations. Every claim below was
checked against the paper (Ashuach et al., *Nature Methods* 2023) and against the implementation in
**scvi-tools 1.4.2** (`scvi/module/_multivae.py`), because the two differ in three places — see §9.

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

Two things to fix in the figure later (§9): "Precision weighted average" is the wrong name for the
fusion, and the arrow from the ATAC decoder into the imputed expression matrix is misleading.

### The letters

| figure | meaning | code name | this document |
|---|---|---|---|
| `Nr`, `Na` | number of RNA cells, number of ATAC cells | `n_obs` | `N` |
| `G` | genes (RNA features) | `n_input_genes` | `G` |
| `P` | peaks (ATAC features) | `n_input_regions` | `P` |
| `B` | latent width — the length of the cell's description | `n_latent` | `L` |
| `H` | the latent matrix, one row per cell | `qz_m` | `z` |

Your figure uses `B` for the latent width. I use `L` below to avoid confusing it with minibatch size,
which is also usually called `B`. Same quantity: 30 in your runs.

---

## 1. The task, in plain words

You have two ways of looking at the same tissue:

- **RNA**: for each cell, how many reads came off each of 36,601 genes.
- **ATAC**: for each cell, how many fragments landed in each of 89,796 open-chromatin regions.

Some cells were assayed both ways. Most were assayed one way. Two things are hard:

1. The two measurements are in **different units** — you cannot put "UMIs of GAPDH" and "fragments in
   peak 41,203" in the same column and compare.
2. Most cells have **only one** of the two, so you cannot simply match them up.

MultiVI solves both at once by learning a **shared short description per cell** — 30 numbers — such that
either measurement can be *predicted back* from those 30 numbers. Cells with only RNA get their 30 numbers
from the RNA; cells with only ATAC get them from the ATAC; cells with both get a blend. And because the
description is enough to predict the other modality, "predict the missing RNA for an ATAC-only cell" becomes
just "run the RNA predictor on that cell's 30 numbers".

That is the whole idea. Everything below is machinery.

---

## 2. The encoder, in plain words

**What goes in.** One cell: a list of 36,601 counts (RNA) or 89,796 counts (ATAC). Typically ~95% are zero.

**What comes out.** Two lists of 30 numbers: a best guess `mu` and an uncertainty `sigma`.

**What it is.** A funnel of two matrix multiplications. The first layer looks at *every* feature at once and
writes 128 intermediate numbers. The second layer mixes those 128 into another 128. Two small final layers
turn those into 30 + 30.

**Where the matrices come from.** Each layer is a table of weights `W` (128 x 36,601 for the first one) plus
128 biases. At the start the numbers are random. They are learned: the model makes a prediction, measures how
wrong it is, and nudges every weight in the direction that reduces the error, a few million weights at a
time. Nothing is hand-designed.

**What it is learning.** Each of the 128 rows of `W1` is a direction in gene space — a pattern of "these
genes go up together". The encoder is learning which patterns matter for describing a cell. Layer two
combines patterns into higher-level ones ("this looks like a B cell"). The 30 outputs are the coordinates of
the cell in that learned space.

**Why also `sigma`.** A cell's 30 numbers are not a fact, they are an estimate. `sigma` says how confident
the estimate is per coordinate. This matters because the model is trained on thousands of cells with very
different depth; a shallow cell should get a vaguer description, and the model needs a way to say so.

**The one-line version.** *The encoder answers: "given this cell's counts, where does it sit in a 30-dimensional
map of cell states, and how sure am I?"*

---

## 3. The sample, in plain words

Between encoder and decoder there is one small step:

> `z = mu + sigma * noise`

Take the guess, jitter it by its own uncertainty. During training the decoder must work from the jittered
position, so it cannot rely on exact coordinates — that is what forces the 30 numbers to form a smooth,
continuous map instead of a lookup table of memorised cells. At inference time the jitter is switched off
and `z = mu`.

---

## 4. The decoder, in plain words

**What goes in.** The 30 numbers.

**What comes out.** 36,601 numbers (RNA) or 89,796 numbers (ATAC) — but **not counts**.

This is the part people misread, so it is worth being blunt: the decoder does not predict "this cell has 4
reads of GAPDH". It predicts **the parameters of a probability distribution** — for each gene, the expected
rate and how over-dispersed it is; for each peak, the probability of being open. The loss then asks "how
likely were the observed counts under that distribution?".

**Why.** Counts are noisy. A gene with a true rate of 3 might read 1, or 5, or 0. Predicting "3" and being
off by 2 is not a mistake, it is sampling. So the well-posed target is the rate, not the realized count.

**The consequence you will see in every result.** A distribution's mean is never exactly zero and never
sparse. So the decoder's output is smooth and dense: applied to your 11,000 ATAC cells it produced a matrix
**22.5x denser** than a measured RNA row. That is not a bug or a threshold; it is what a conditional mean
looks like. It is also why the GRN step behaves the way it does (§8).

**The one-line version.** *The decoder answers: "given a position on the map, what count distribution should
each gene / each peak have?"*

---

## 5. Why there are two of each

- **Two encoders** because RNA and ATAC are different data types with different noise: RNA is
  over-dispersed counts, ATAC is essentially binary. One encoder cannot serve both well.
- **Two decoders** for the same reason, in reverse.
- **Masking** is what lets one model hold both: a cell with no ATAC counts contributes zero ATAC
  reconstruction error. The mask is computed, not supplied: `mask = (rowsum of that modality > 0)`.

So a cell with only ATAC gets its `z` from the ATAC encoder alone, and is only ever asked to explain its
peaks — never its (missing) genes. Later, the RNA decoder can still be run on that `z` to guess the genes.
That single asymmetry is the imputation.

---

## 6. The mathematics

Notation: modality `m` in `{R, A}`; cell `i`; `x_i^R` its gene counts (length `G`), `x_i^A` its peak counts
(length `P`); `L = 30`; `H = 128`.

### 6.1 Encoder (per modality, separate weights)

For input `x` (a minibatch, `B x G` or `B x P`):

```
h1 = LayerNorm( LeakyReLU( W1 x + b1 ) )        W1: H x G  (or H x P)       B x H
h2 = LayerNorm( LeakyReLU( W2 h1 + b2 ) )       W2: H x H                    B x H
mu    = W_mu h2 + b_mu                          W_mu: L x H                  B x L
sigma2 = exp( W_v h2 + b_v )                    W_v:  L x H                  B x L
```

Exact details from the code (`scvi.nn.Encoder`, used at `_multivae.py:367`):

- activation is **LeakyReLU**, not ReLU;
- **LayerNorm** after each layer (`use_layer_norm="both"`);
- dropout 0.1 inside the layers;
- the variance is `exp(linear)`, so it is positive by construction (`var_activation=torch.exp`,
  `var_eps=0`) — the code stores a **variance**, not a standard deviation.

There is a second, parallel encoder for library size:

```
l = softplus( w_l h2' + b )                     length 1 per cell
```

It answers "how deep was this cell sequenced?", so the decoder can separate "more reads" from "more biology".

### 6.2 Sampling (reparameterisation)

```
eps ~ Normal(0, I)                B x L
z   = mu + sqrt(sigma2) * eps     B x L
```

Written this way, the randomness is outside the network, so gradients flow through `mu` and `sigma` normally.

### 6.3 Fusion — one latent from two beliefs

For a cell with both modalities you now have two guesses. They are combined by a **masked weighted mean**:

```
W_m = masked_softmax(w_m , mask_m)              weights across modalities, sum to 1, absent -> 0
mu  = sum_m  W_m * mu_m                         B x L
var = sum_m  sqrt(W_m) * var_m                  B x L
```

where `w` is one number per modality (`equal` = 0.5/0.5, the default), or a learned global pair
(`universal`), or a learned per-cell table (`cell`).

Three things to notice:

- A cell with one modality keeps that encoder's belief exactly: `W` is 1 for the present modality, 0 for the
  absent one.
- This is **not** precision weighting and **not** a product of experts. Nothing divides by a variance. The
  paper describes it the same way: `Z_c = w_A Z_c^A + w_R Z_c^R`, `sum_m w_m = 1`.
- The variance line is a heuristic pool (`sqrt` of the weights), exact only when a cell has one modality.

### 6.4 RNA decoder

```
h  = LeakyReLU stack( z )                       B x H   (two layers, LayerNorm)
px_scale = softmax( W_s h + b_s )               B x G   proportions over genes, sum to 1 per cell
rate     = exp(l) * px_scale                    B x G   expected counts
theta    = exp( W_t h + b_t )                   G       per-gene inverse dispersion
pi       = W_d h + b_d                          B x G   dropout logits
```

Then the likelihood, per cell, summed over genes:

```
ZINB (default):
  -log p(x | rate, theta, pi) = -sum_g log[ ZINB(x_g ; rate_g, theta_g, pi_g) ]

NB (if gene_likelihood="nb"): drop pi
Poisson: rate only
```

`softmax` is why `px_scale` is a *within-cell proportion*: it says what fraction of that cell's reads went to
each gene, and multiplying by the cell's library size `exp(l)` turns it into expected counts. That is also
the number your pipeline reports as the imputed expression, after `log1p(1e4 * px_scale)`.

### 6.5 ATAC decoder

```
h = LeakyReLU stack( z )                        B x H
p = sigmoid( W_p h + b_p )                      B x P   probability each peak is open
r = sigmoid(region_factors)                     P       per-peak accessibility prior, learned
```

Loss per cell, summed over peaks, against the binarised counts:

```
BCE = -sum_j [ x_j log(p_j * l * r_j) + (1 - x_j) log(1 - p_j * l * r_j) ] ,  x_j = 1[count > 0]
```

### 6.6 The objective

Per minibatch, averaged over cells:

```
L = mean_i [  mask_R(i) * (-log p(x_i^R | z_i))          reconstruction, RNA
            + mask_A(i) * BCE(x_i^A , p_i)               reconstruction, ATAC
            + beta * KL( q(z_i | x_i) || Normal(0, I) )  keep the map tidy
            + Jeffreys( q_R(z_i) || q_A(z_i) ) ]         make the two encoders agree
```

with `beta = 1` and no separate weight on the Jeffreys term in scvi-tools 1.4.2 — it is added to the KL
directly (`_multivae.py:868`).

The **KL** term, written out:

```
KL = 0.5 * sum_k [ mu_k^2 + sigma_k^2 - log(sigma_k^2) - 1 ]
```

It pulls every cell's 30 numbers toward a standard normal: centred at zero, unit spread, no clumps. Without
it the model could give every cell its own private coordinate and the map would be meaningless.

The **Jeffreys** term (symmetric KL between the two modality beliefs), for cells with both:

```
J = KL(q_R || q_A) + KL(q_A || q_R)
  = 0.5 * sum_k [ sigma_R^2/sigma_A^2 + sigma_A^2/sigma_R^2 - 2
                  + (mu_R - mu_A)^2 * (1/sigma_R^2 + 1/sigma_A^2) ]
```

Read it as a tug-of-war: the means must coincide, and the uncertainties must match. It is computed only where
a cell has both modalities (`mask_R & mask_A`, `_multivae.py:934-970`). **With no co-assayed cells this term
is identically zero**, and nothing else forces the two encoders into the same space — which is exactly the
failure mode measured earlier (FOSCTTM 0.43 when paired rows were presented separately, 0.003 when collapsed).

### 6.7 How it is optimised

Shuffled minibatches of 128 cells, Adam, 200 epochs, early stopping on a held-out ELBO. One pass over 11,000
cells is 86 gradient steps. Nothing is solved in closed form; every matrix above is updated by gradient
descent.

---

## 7. Where every matrix is, and how big

At your widths: `G = 36,601`, `P = 89,796`, `H = 128`, `L = 30`. Counts measured by instantiating the module.

| module | layers | parameters |
|---|---|---|
| `z_encoder_expression` | 36,601 -> 128, 128 -> 128, 128 -> 30 (mean), 128 -> 30 (var) | 4,709,308 |
| `l_encoder_expression` | 36,601 -> 128, 128 -> 128, 128 -> 1 | 4,701,697 |
| `z_decoder_expression` | 30 -> 128, 128 -> 128, then 3 heads 128 -> 36,601 (scale, dispersion, dropout) | 14,185,067 |
| `z_encoder_accessibility` | 89,796 -> 128, 128 -> 128, 128 -> 30 x 2 | 11,518,268 |
| `l_encoder_accessibility` | 89,796 -> 128, 128 -> 128, 128 -> 1 | 11,510,657 |
| `z_decoder_accessibility` | 30 -> 128, 128 -> 128, 128 -> 89,796 (sigmoid) | 11,604,164 |
| protein stubs (unused, 0 proteins) | — | 85,948 |
| **total** | | **58,441,506** |

Two facts that follow directly from the table:

- **The first layer dominates.** `W1` for the ATAC encoder alone is 11.5M numbers — 20% of the entire model in
  one matrix, because it must map 89,796 peaks to 128 numbers.
- **RNA pays for three output heads** (mean, dispersion, dropout) while ATAC pays for one. That is why 36,601
  genes cost 14.2M in the decoder but 89,796 peaks only cost 11.6M.

---

## 8. One cell, end to end, with numbers

Take one ATAC-only cell from your 11,000.

**Encode.** Its row of 89,796 peak counts (about 2,500 non-zero) goes through `W1` (128 x 89,796): ~11.5M
multiply-adds. Two more small layers. Out: `mu` (30 numbers), `sigma2` (30 numbers). The RNA encoder is
skipped in effect — the cell has no RNA, so its mask is 0 and it contributes nothing there.

**Fuse.** One belief, so `W_R = 0`, `W_A = 1`. `z = mu`.

**Decode.** `z` (30 numbers) goes into the RNA decoder: 30 -> 128 -> 128 -> 36,601, three times (scale,
dispersion, dropout): ~14M multiply-adds. Out: `px_scale`, a 36,601-long vector of proportions summing to 1.

**That vector is the imputed expression for that cell.** Your pipeline applies `log1p(1e4 * px_scale)` and
stores it as one column of the `36,601 x 11,000` imputed matrix, then stacks it under the measured RNA rows to
make the `(Nr + Na) x G` matrix that Arboreto reads.

**Why it is smooth.** Every entry is a conditional mean from a softmax. It has no zeros, and neighbouring
cells (nearby `z`) get near-identical rows. Measured: 22.5x denser than a real row, gene means off by 0.336 sd.

**Why the GRN result follows.** GRNBoost2 looks for co-variation across cells. A matrix whose imputed half is
a smooth function of a 30-dimensional latent has strong, smooth co-variation that is *not* in the data — the
trees will happily split on it. That is the mechanism behind AUPRC 0.006 against prevalence 0.005, and early
precision 0.000: not a broken network code, but an input whose imputed half is too regular.

---

## 9. Three places the figure and the code disagree

1. **"Posterior fusion: Precision weighted average".** The code does a masked weighted mean of the means
   (`masked_softmax` then weighted sum; variances weighted by `sqrt(w)`). There is no `1/sigma^2` anywhere.
   The paper agrees with the code: `Z_c = w_A Z_c^A + w_R Z_c^R`. Relabel the box **"Posterior fusion —
   masked weighted mean"**.
2. **The arrow from the ATAC decoder into the imputed expression matrix.** The imputed expression matrix
   (`G x Na`) comes from the **RNA decoder** applied to ATAC cells' `z`. The ATAC decoder reconstructs peaks;
   its output is a `P x Na` peak matrix that your pipeline never uses. Draw that arrow out of the picture, or
   route it to a separate "imputed peaks (not used)" box.
3. **`n_hidden` default.** It is `min(128, sqrt(n_features))`, so it is 128 for any realistic dataset — it does
   not grow with the peak count. The 469.8M-parameter scenario only happens if you *set* `h=500`, not if you
   leave it unset.

One more that is not in the figure but matters for reading the results: the **default likelihood is ZINB**, not
NB, so the dropout head (`pi`) is actually used. If you ever set `gene_likelihood="nb"`, that head is allocated
and ignored — a 4.7M-parameter no-op.

---

## 10. Glossary

| term | plain meaning |
|---|---|
| encoder | funnel from a cell's counts to 30 numbers + 30 uncertainties |
| decoder | reverse funnel from 30 numbers to distribution parameters for every feature |
| latent `z` | the 30 numbers; the cell's coordinates on the learned map |
| `mu`, `sigma2` | per-coordinate guess and variance from one encoder |
| mask | `rowsum(counts) > 0`; which modality a cell actually has |
| fusion | combining the RNA and ATAC beliefs into one `z` |
| `px_scale` | per-cell proportions over genes (softmax output); the imputed expression |
| `rate` | `px_scale` x library size; expected counts fed to the likelihood |
| dispersion `theta` | how over-dispersed a gene's counts are; NB/ZINB shape |
| dropout `pi` | ZINB zero-inflation logit; probability a zero is technical |
| `p` | per-peak probability of being open (ATAC decoder output) |
| region factor | learned per-peak scaling; how accessible that peak is overall |
| KL | penalty keeping the map centred and smooth |
| Jeffreys | penalty forcing the two encoders to agree on co-assayed cells |
| ELBO | the objective being maximised = likelihood of the data minus the penalties |
