#!/usr/bin/env python3
import yaml, scanpy as sc

cfg = yaml.safe_load(open("config.yaml"))
adata = sc.read_h5ad(cfg["data_path"])
sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[adata.obs.n_genes_by_counts > cfg["min_genes_per_cell"], :]
adata = adata[adata.obs.pct_counts_mt < cfg["max_mito_pct"], :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write("data/processed/01_preprocessed.h5ad")
print("✅ Preprocessing complete.")
