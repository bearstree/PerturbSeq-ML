# PerturbSeq-ML

Single-cell CRISPRi Perturb-seq analysis of public K562 data, from quality control and transcriptional-state analysis to perturbation responses, classification, and AnnData export.

## Data and workflow

Source: Replogle et al. (2022), *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq*. [Public data](https://doi.org/10.25452/figshare.plus.20029387) · [Paper](https://doi.org/10.1016/j.cell.2022.05.013). Download `K562_essential_raw_singlecell_01.h5ad` to `data/raw/`. The ~10 GB raw file is excluded from Git. The original uses a legacy AnnData encoding; `scripts/prepare_k562_subset.py` creates a compatible sparse analysis H5AD without changing it.

Run the notebooks in order from the repository root:

| Notebook | Analysis |
| --- | --- |
| [01](notebooks/01_data_import_QC.Rmd) | Import, sparse counts, full-data QC, seeded coverage subset |
| [02](notebooks/02_standard_scRNA_analysis.Rmd) | HVGs, PCA, UMAP, clustering, K562 states |
| [03](notebooks/03_perturbation_response.Rmd) | Powered perturbation cohort, target expression, exploratory DE, aggregate effects, E-distance |
| [04](notebooks/04_perturbation_embeddings.Rmd) | Cosine similarity, perturbation PCA, descriptive response clusters |
| [05](notebooks/05_ML_perturbation_classifier.Rmd) | Held-out-cell classification of known classes |
| [06](notebooks/06_AnnData_Scanpy_interop.Rmd) | H5AD export and R/Python round-trip |

Rendered HTML for each notebook is in `notebooks/` (Notebook 01 uses `.nb.html`).

## Setup

Use R with Seurat, anndataR, Matrix, dplyr, ggplot2, pheatmap, ggrepel, scperturbR, reticulate, and rmarkdown. Python needs anndata, h5py, scipy, numpy, pandas, scikit-learn, and joblib; Scanpy is optional. These notebooks were executed with R 4.6.1 and Python 3.14 on Windows. Notebook 03 and 05–06 point reticulate at `C:/Python314/python.exe`; change that path for another environment.

```powershell
python scripts/prepare_k562_subset.py
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' -e "for (n in sprintf('notebooks/%02d_%s.Rmd', 1:6, c('data_import_QC','standard_scRNA_analysis','perturbation_response','perturbation_embeddings','ML_perturbation_classifier','AnnData_Scanpy_interop'))) rmarkdown::render(n, knit_root_dir=getwd())"
```

## Executed results

- Full source: 310,385 cells × 8,563 genes. QC retained 305,906 cells. The seeded 12,000-cell coverage subset retained all 2,058 perturbation labels; clustering found 10 descriptive K562 states.
- The focused Part 2 cohort has 12 adequately represented targets (245–249 cells each) and 493 non-targeting controls. It supports detailed comparisons that the broad coverage subset cannot.
- Exploratory cell-level DE found 94 DONSON and 2 BUB1 hits under the notebook's adjusted-p-value and effect-size cutoffs. GAB2 had the largest aggregate effect norm (24.474) and E-distance from controls (6.484). DONSON–RRM2 had the highest effect-vector cosine similarity (0.696).
- The 13-class multinomial elastic-net classifier used 2,766 training and 692 test cells. Held-out accuracy was **0.410**, macro-F1 **0.386**, versus a **0.143** majority-class accuracy baseline. This is moderate, uneven performance despite exceeding baseline.
- The focused Seurat object exported to H5AD and passed R and Python AnnData checks (3,458 cells × 8,563 genes, 13 labels, PCA embedding).

## Interpretation and limits

K562 clusters are transcriptional/cellular states, not immune-cell types. The 12,000-cell subset measures label coverage; it is too sparse per target for the detailed Part 2 analyses. The focused cohort was sampled separately from the public raw file with seed 42 and the Part 1 QC limits.

Individual cells are not independent biological replicates, so cell-level DE is exploratory. One aggregate per label summarizes expression but does not provide biological-replicate-aware inference. E-distance and E-tests with 199 cell-label permutations are exploratory. Detailed perturbation analysis covers 12 targets rather than all 2,058 labels. The classifier's random cell split tests new cells from **known** classes; it does not test unseen targets or guide-wise generalization. Most targeting labels have one observed guide pair, and substantial class overlap remains. Results do not establish experimental validation or a causal mechanism.

## Reproducibility notes

Scanpy was unavailable in the execution environment, so Scanpy-specific compatibility was not directly tested. R and Python AnnData round-trip checks passed. Seurat sanitized underscores to hyphens in two duplicate gene names; raw counts remained aligned. The focused object has PCA but no UMAP embedding. `artifacts/manifest.csv` records selected output paths, sizes, and hashes; `artifacts/sessionInfo.txt` records the R environment.

## Repository files

`notebooks/` contains R Markdown sources and executed HTML; `results/` contains lightweight tables; `models/` contains the fitted classifier; `artifacts/` contains the export, manifest, and R session information. `data/raw/` and large `data/processed/` files remain local and are excluded from Git. The export H5AD is also local only; its dimensions and metadata checks are visible in Notebook 06.

## License

Repository code and documentation are available under the [MIT License](LICENSE). The public dataset is provided by its original source and is not relicensed by this repository.
