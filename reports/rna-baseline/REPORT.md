# RNA-only baseline — GRNBoost2 on stacked PBMC 3k + PBMC 4k RNA

**Baseline:** no integration, no imputation — two RNA gene-expression matrices stacked and handed directly to Arboreto.
**GRN:** Arboreto GRNBoost2 (TRRUST TFs as regulators, TRRUST-only targets)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Inputs

- **3k RNA:** PBMC multiome 3k granulocyte-sorted, Gene Expression (36601 genes, GRCh38-2020-A), 2711 cells.
- **4k RNA:** 10x PBMC 4k from a healthy donor (v2 / Cell Ranger 2.1.0, 33694 genes, GRCh38-3.0.0), 4340 cells.
- **Stacking:** gene symbols common to both references (**21932** genes); duplicate symbols summed per reference. Each matrix log1p(CPM)-normalized independently, then cells stacked (3k rows first, then 4k).
- **All-cells matrix:** 7051 cells x 21932 genes.

## 2. Arboreto GRNBoost2 (TRRUST regulators, TRRUST-only targets)

- Cells: 7051
- Regulators (TRRUST TFs present in expression): 2808
- Target genes (TRRUST TFs present): 2808
- Edges inferred: 751,233

**Top 25 edges:**

| TF -> target | importance |
|---|---|
| CD8A -> CD8B | 258.4 |
| CD8B -> CD8A | 203.0 |
| HLA-B -> HLA-C | 116.3 |
| HLA-C -> HLA-B | 112.0 |
| FOS -> DUSP1 | 109.5 |
| S100A6 -> S100A4 | 104.7 |
| CD79A -> MS4A1 | 100.5 |
| HLA-B -> HLA-A | 99.8 |
| TPT1 -> RPS3A | 97.9 |
| DUSP1 -> FOS | 97.6 |
| MS4A1 -> CD79A | 96.4 |
| RPS3A -> TPT1 | 92.7 |
| LYZ -> S100A9 | 90.8 |
| GNLY -> PRF1 | 90.2 |
| HLA-B -> B2M | 84.8 |
| FOS -> JUN | 84.6 |
| S100A4 -> S100A6 | 82.3 |
| PRF1 -> GNLY | 81.7 |
| RPS3A -> ALDOA | 79.9 |
| HLA-DRB1 -> HLA-DQB1 | 79.9 |
| HLA-DRA -> HLA-DRB1 | 78.2 |
| RPS3A -> RPL3 | 76.5 |
| BCL2 -> CDK6 | 75.4 |
| MNDA -> S100A9 | 73.9 |
| LYZ -> VCAN | 73.8 |

## 3. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=1433 (1433/8751, 16.4%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0100 | 0.0001 |
| 500 | 0.0040 | 0.0002 |
| 1000 | 0.0050 | 0.0006 |
| 5000 | 0.0054 | 0.0031 |
| 10000 | 0.0042 | 0.0048 |
| 751233 | 0.0019 | 0.1638 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=1841 (1841/96846, 1.9%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0100 | 0.0000 |
| 500 | 0.0040 | 0.0000 |
| 1000 | 0.0030 | 0.0000 |
| 5000 | 0.0038 | 0.0002 |
| 10000 | 0.0036 | 0.0004 |
| 751233 | 0.0025 | 0.0190 |

---
