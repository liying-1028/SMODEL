c# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:57:23 2024

@author: liying
"""

import pandas as pd
import scanpy as sc
import anndata
import sklearn
import numpy as np
import scanpy as sc
import pandas as pd
from typing import Optional


def clr_normalize_each_cell(adata, inplace=True):
    """Normalize count vector for each cell, i.e. for each row of .X"""

    import numpy as np
    import scipy

    def seurat_clr(x):
        # TODO: support sparseness
        s = np.sum(np.log1p(x[x > 0]))
        exp = np.exp(s / len(x))
        return np.log1p(x / exp)

    if not inplace:
        adata = adata.copy()

    # apply to dense or sparse matrix, along axis. returns dense matrix
    adata.X = np.apply_along_axis(
        seurat_clr, 1, (adata.X.A if scipy.sparse.issparse(adata.X) else np.array(adata.X))
    )
    return adata

def lsi(
        adata: anndata.AnnData, n_components: int = 20,
        use_highly_variable: Optional[bool] = None, **kwargs
) -> None:
    r"""
    LSI analysis (following the Seurat v3 approach)
    """
    if use_highly_variable is None:
        use_highly_variable = "highly_variable" in adata.var
    adata_use = adata[:, adata.var["highly_variable"]] if use_highly_variable else adata
    X = tfidf(adata_use.X)
    # X = adata_use.X
    X_norm = sklearn.preprocessing.Normalizer(norm="l1").fit_transform(X)
    X_norm = np.log1p(X_norm * 1e4)
    X_lsi = sklearn.utils.extmath.randomized_svd(X_norm, n_components, **kwargs)[0]
    X_lsi -= X_lsi.mean(axis=1, keepdims=True)
    X_lsi /= X_lsi.std(axis=1, ddof=1, keepdims=True)
    # adata.obsm["X_lsi"] = X_lsi
    adata.obsm["X_lsi"] = X_lsi[:, 1:]
    
def tfidf(X):
    r"""
    TF-IDF normalization (following the Seurat v3 approach)
    """
    idf = X.shape[0] / X.sum(axis=0)
    if scipy.sparse.issparse(X):
        tf = X.multiply(1 / X.sum(axis=1))
        return tf.multiply(idf)
    else:
        tf = X / X.sum(axis=1, keepdims=True)
        return tf * idf

# read data
file_fold = '/data/liying/work/data/' #please replace 'file_fold' with the download path

adata_omics1 = sc.read_h5ad(file_fold + 'adata_RNA.h5ad')
adata_omics2 = sc.read_h5ad(file_fold + 'adata_ADT.h5ad')
adata_omics3 = sc.read_h5ad(file_fold + 'adata_ATAC.h5ad')

adata_omics1.var_names_make_unique()
adata_omics2.var_names_make_unique()
adata_omics3.var_names_make_unique()


n_protein = adata_omics2.n_vars

sc.pp.highly_variable_genes(adata_omics1, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata_omics1, target_sum=1e4)
sc.pp.log1p(adata_omics1)

adata_omics1_high =  adata_omics1[:, adata_omics1.var['highly_variable']]

gene3000 = adata_omics1_high.to_df()
gene3000.to_csv('/data/liying/work/data/results/RNA.csv', sep=',', index=True, header=True)


# Protein
adata_omics2 = clr_normalize_each_cell(adata_omics2)
mmtv_ADT = adata_omics2.to_df()
mmtv_ADT.to_csv('/data/liying/work/data/results/ADT.csv', sep=',', index=True, header=True)


# ATAC
lsi(adata_omics3, use_highly_variable=False, n_components=n_protein + 1)#revise
adata_omics3.obsm['feat'] = adata_omics3.obsm['X_lsi'].copy()
feat_data = adata_omics3.obsm['feat']
feat_df = pd.DataFrame(feat_data, index=adata_omics3.obs_names)
output_file = "/data/liying/work/data/results/ATAC_lsi.csv"
feat_df.to_csv(output_file, sep=',', index=True, header=True)

