"""Full-dataset QC and seeded, perturbation-stratified H5AD conversion."""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/raw/K562_essential_raw_singlecell_01.h5ad")
    parser.add_argument("--output", default="data/processed/01_k562_analysis.h5ad")
    parser.add_argument("--max-cells", type=int, default=12000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    source = ad.read_h5ad(args.input, backed="r")
    obs = source.obs.copy()
    names = source.var["gene_name"].astype(str).to_numpy()
    mt = np.char.startswith(names.astype(str), "MT-")
    n_count = np.zeros(source.n_obs)
    n_feature = np.zeros(source.n_obs, dtype=int)
    mt_count = np.zeros(source.n_obs)
    for start in range(0, source.n_obs, 1024):
        stop = min(start + 1024, source.n_obs)
        x = np.asarray(source.X[start:stop, :])
        n_count[start:stop] = x.sum(axis=1)
        n_feature[start:stop] = np.count_nonzero(x, axis=1)
        if mt.any():
            mt_count[start:stop] = x[:, mt].sum(axis=1)
        if start % 32768 == 0:
            print(f"QC: {stop}/{source.n_obs}", flush=True)
    obs["nCount_RNA"] = n_count
    obs["nFeature_RNA"] = n_feature
    if mt.any():
        obs["percent.mt"] = np.divide(100 * mt_count, n_count, out=np.zeros_like(n_count), where=n_count > 0)
    fb = np.quantile(n_feature, [0.005, 0.995])
    cb = np.quantile(n_count, [0.005, 0.995])
    keep = (n_feature >= fb[0]) & (n_feature <= fb[1]) & (n_count >= cb[0]) & (n_count <= cb[1])
    if mt.any():
        keep &= obs["percent.mt"].to_numpy() <= 20
    rng = np.random.default_rng(args.seed)
    kept = np.flatnonzero(keep)
    labels = obs.iloc[kept]["gene"].astype(str).to_numpy()
    groups = {label: kept[labels == label] for label in np.unique(labels)}
    if len(kept) > args.max_cells:
        if len(groups) > args.max_cells:
            raise ValueError("max-cells cannot preserve every perturbation")
        per_target = min(5, args.max_cells // len(groups))
        guaranteed = np.concatenate([rng.choice(v, min(len(v), per_target), replace=False) for v in groups.values()])
        remaining = np.setdiff1d(kept, guaranteed)
        extra = rng.choice(remaining, args.max_cells - len(guaranteed), replace=False)
        selected = np.sort(np.concatenate([guaranteed, extra]))
    else:
        selected = kept
        per_target = 0
    chunks = []
    for batch in np.array_split(selected, max(1, int(np.ceil(len(selected) / 512)))):
        chunks.append(sparse.csr_matrix(np.asarray(source.X[batch, :])))
    var = source.var.copy()
    # Preserve unique gene symbols while retaining Ensembl IDs in var.
    var["ensembl_id"] = var.index.astype(str)
    seen = {}
    symbols = []
    for name in names:
        count = seen.get(name, 0)
        symbols.append(name if count == 0 else f"{name}_{count}")
        seen[name] = count + 1
    var.index = symbols
    result = ad.AnnData(X=sparse.vstack(chunks, format="csr"), obs=obs.iloc[selected].copy(), var=var)
    result.uns["qc"] = {"before": int(source.n_obs), "after": int(keep.sum()), "analysis": int(len(selected)), "feature_lower": float(fb[0]), "feature_upper": float(fb[1]), "count_lower": float(cb[0]), "count_upper": float(cb[1]), "mt_gene_count": int(mt.sum()), "mt_max_percent": 20.0 if mt.any() else -1.0, "seed": args.seed, "guaranteed_per_target": per_target}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(args.output, compression="gzip")
    Path("results").mkdir(exist_ok=True)
    pd.DataFrame([result.uns["qc"]]).to_csv("results/01_qc_report.csv", index=False)
    obs.loc[keep, "gene"].value_counts().rename_axis("gene").reset_index(name="n_cells").to_csv("results/01_perturbation_cell_counts.csv", index=False)
    result.obs["gene"].value_counts().rename_axis("gene").reset_index(name="n_cells").to_csv("results/01_analysis_perturbation_counts.csv", index=False)
    print(result.uns["qc"])
    print("Sparse matrix:", result.X.shape, result.X.nnz)
    source.file.close()


if __name__ == "__main__":
    main()
