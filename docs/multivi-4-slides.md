# MultiVI in 4 slides

Speaker notes only — one line per bullet is the whole slide. Numbers are the measured ones from
`docs/multivi-architecture-and-scaling.md` (repo width: 36,601 genes, 89,796 peaks, h=128, L=30).

---

## Slide 1 — What MultiVI is

- Two encoders, **one shared latent**: RNA -> `z`, ATAC -> `z`, fused by a precision-weighted average.
- Two decoders, each reconstructing only where data exist — the modality mask is `rowsum > 0`, not user input.
- Paired cells get a **Jeffreys penalty** that forces the two encoders to agree. No paired cells, no penalty.
- Cost is `5hG + 3hP` parameters — **58.4 M** here, 0.94 GB with Adam.
- **There is no `N x N` object anywhere in the model.** That single fact decides everything on slide 2.

*Figure:* the nine boxes below, left to right, with the tensor on every arrow.

### Slide 1 boxes

Nine boxes in four stages. Two parallel columns on the left (RNA top, ATAC bottom), one trunk in the
middle, two parallel columns on the right, loss underneath.

```
        RNA                                        ATAC
  ┌───────────────┐                          ┌────────────────┐
  │ 1 RNA counts  │ N x G sparse             │ 2 ATAC counts  │ N x P sparse
  └───────┬───────┘                          └───────┬────────┘
          │ x_r  B x G                               │ x_a  B x P
  ┌───────▼───────┐                          ┌───────▼────────┐
  │ 3 RNA encoder │ G -> h -> h              │ 4 ATAC encoder │ P -> h -> h
  └───────┬───────┘                          └───────┬────────┘
          │ mu_r, sigma_r  B x L                     │ mu_a, sigma_a  B x L
          │            ╲        (Jeffreys)  ╱        │
          │             ╲   only on paired ╱         │
          │              ▼                 ▼         │
          │        ┌──────────────────────────┐      │
          └───────▶│ 5 Posterior fusion       │◀─────┘
                   │   precision-weighted avg │   + per-cell weights w
                   └────────────┬─────────────┘
                                │ mu, sigma  B x L
                   ┌────────────▼─────────────┐
                   │ 6 Sample z               │
                   │   z = mu + sigma * eps   │
                   └────────────┬─────────────┘
                                │ z  B x L        --> also the integration output, N x L
              ┌─────────────────┴──────────────────┐
     ┌────────▼────────┐                  ┌────────▼────────┐
     │ 7 RNA decoder   │ h -> 3 heads     │ 8 ATAC decoder  │ h -> Bernoulli
     └────────┬────────┘                  └────────┬────────┘
              │ px_scale  B x G                    │ p  B x P
     ┌────────▼────────────────────────────────────▼────────┐
     │ 9 Loss    RNA NB recon + ATAC BCE recon + KL + Jeffreys│
     └──────────────────────────────────────────────────────┘
```

| # | box label (write this) | second line inside the box | what enters |
|---|---|---|---|
| 1 | RNA counts | N x G sparse, 36,601 genes | — |
| 2 | ATAC counts | N x P sparse, 89,796 peaks | — |
| 3 | RNA encoder | G->h->h, mu_r / sigma_r : B x L | x_r |
| 4 | ATAC encoder | P->h->h, mu_a / sigma_a : B x L | x_a |
| 5 | Posterior fusion | precision-weighted mean, mu / sigma : B x L | mu_r, sigma_r, mu_a, sigma_a, w |
| 6 | Sample z | z = mu + sigma * eps, B x L | mu, sigma |
| 7 | RNA decoder | h->G, three heads (mean, dispersion, dropout) | z |
| 8 | ATAC decoder | h->P, Bernoulli probability | z |
| 9 | Loss | recon RNA + recon ATAC + KL + Jeffreys | px_scale, p, mu, sigma |

Three dashed callouts, not boxes:

- **masks** — from boxes 1 and 2 straight to box 9, labelled `m = (rowsum > 0)`: a cell with no counts in
  a modality contributes no reconstruction error there.
- **Jeffreys** — a double arrow between boxes 3 and 4, labelled *paired cells only*. If you draw one
  thing on this slide, draw this: it is the only term that ties the two modalities together.
- **w** — into box 5, labelled `equal | universal | per-cell`.

Two outputs, drawn as plain arrows off the diagram, not boxes:

- `z` off box 6 → the joint embedding, `N x L`, what the integration metrics score.
- `px_scale` off box 7 → the imputed expression matrix, `n_atac x G`, what feeds GRNBoost2.

If nine boxes is too many for one slide, collapse to six: `RNA counts | ATAC counts` → `Encoders` (boxes
3+4) → `Posterior fusion` → `Sample z` → `Decoders` (boxes 7+8) → `Loss`. The Jeffreys callout must
survive the collapse; it is the one thing a reader cannot guess.

---

## Slide 2 — Does it scale to 10^6 cells?

| | now (11k cells) | 1M cells |
|---|---|---|
| per-cell cost | 1.6 ms/epoch (measured, 8 CPU threads) | same |
| one 200-epoch fit | 58 min | 44 h CPU / ~1 h GPU |
| model + Adam | 0.94 GB | 1.9 GB |
| decoded matrix | 1.6 GB | **146 GB** |
| `N x N` terms | none | none |

- Bottleneck moves with scale: at 10^4-10^5 it is **compute** (`P` is the knob); at 10^6 it is the
  **sparse-to-dense data path** and the **export**, not the model.
- Peaks and genes are not symmetric: ATAC costs `3hP` (34.6 M at 89,796 peaks), RNA costs `5hG`
  (23.6 M at 36,601 genes) — the peak axis is 2.5x wider but only 1.5x the parameters, because RNA
  pays for three decoder heads.
- `n_hidden` defaults to `sqrt(n_regions)` — leave it unset at 250k peaks and the model is 3.9x bigger.
- **No published MultiVI result exceeds 89,655 cells.** The ">1 M" claim is the scvi-tools doc page.

---

## Slide 3 — MultiVI vs scSAGA

| | scSAGA | MultiVI |
|---|---|---|
| input | 50-d PCA | raw counts, 126k features |
| training | none, one shot | 200 epochs x SGD |
| time | O(S·N·K), S=25 | O(E·N·h·(G+P)) |
| memory | O(N·K) dense coupling | O(nnz) |
| measured, 11k cells | ~3 min | ~58 min (one 200-epoch fit) |
| 1M cells | 4 TB coupling → impossible | ~44 h CPU / ~1 h GPU |

- Algorithmically it is **cells² vs features x epochs**. scSAGA never touches the peak matrix; MultiVI touches every peak every epoch.
- scSAGA is ~20-40x cheaper at 10^4 — but its work is graph traversal (latency-bound), MultiVI's is dense GEMM (~10^11 FLOP/s). A GPU accelerates one and not the other, which is why the ranking flips.
- scSAGA needs **no paired cells at all** (Gromov-Wasserstein preserves within-domain distances). MultiVI's prior is "paired cells agree" — weaker, but linear.

---

## Slide 4 — What the 11k run actually showed

- **Integration: MultiVI wins.** FOSCTTM 0.008 vs 0.241 (random 0.5). Its two encoders agree; scSAGA's two rows do not.
- **Imputation: both fail the mean control.** MultiVI MSE 0.172, scSAGA 0.124, mean predictor 0.097. MultiVI is 22.5x over-dense and shifts gene means by 0.336 sd; scSAGA preserves means (0.013 sd) but its MSE *equals* the control — no cell-specific signal.
- **GRN: at the noise floor.** AUPRC 0.006 / 0.010 against prevalence 0.005 / 0.006. Hub overlap 0.000 for both. scSAGA wins every precision metric; MultiVI is only more seed-stable (0.546 vs 0.405).
- **The pattern:** better latent ≠ better network. The step that decides the GRN is imputation density, and scSAGA's cruder, marginal-preserving average transfers to GRNBoost2 better than MultiVI's smooth decode.
- **Blocker:** 12 metrics are `n/a` for want of `labels.txt`. Nothing biological is verified yet.

*Next:* add labels, then the 4-config ablation in `report_extrapolated.md` §7.
