# https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
import scanpy as sc
import anndata as ad
import sys
import hdf5plugin
import json

if __name__=='__main__':
    file = sys.argv[1]

    groups_file_path = file.replace("data.h5", "groups.json")
    
    groups_data = []
    with open(groups_file_path) as f:
        groups_data = json.load(f)

    adatas = {}
    for group in groups_data:
        sample_adata = sc.read_10x_h5(group["files"][0])
        sample_adata.var_names_make_unique()
        adatas[group["name"]] = sample_adata

    # Load the data
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()
    print(adata.obs["sample"].value_counts())
    print(adata)

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
    output_folder = file.replace("data.h5", "scanpy"),
    sc.settings.figdir = output_folder
    

    # VIOLIN PLOT
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=".png"
    )
    
    # SCATTER PLOT
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save=".png")
    
    # Save data used for those two plot to a JSON file
    violin_data = {
        'n_genes_by_counts': adata.obs['n_genes_by_counts'].tolist(),
        'total_counts': adata.obs['total_counts'].tolist(),
        'pct_counts_mt': adata.obs['pct_counts_mt'].tolist(),
        'type': 'violin'
    }
    with open(output_folder+'/violin_and_scatter_plot_data.json', 'w') as f:
        json.dump(violin_data, f)

    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.scrublet(adata, batch_key="sample")

    # Saving count data
    adata.layers["counts"] = adata.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)

    # HIGHLY VARIABLE GENES
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    sc.pl.highly_variable_genes(adata, save=".png")
    highly_variable_genes_data = {
        'highly_variable': adata.var['highly_variable'].tolist(),
        'means': adata.var['means'].tolist(),
        'dispersions': adata.var['dispersions'].tolist(),
        # 'dispersion_norm': adata.var['dispersion_norm'].tolist()
    }
    with open(output_folder+'/highly_variable_genes_data.json', 'w') as f:
        json.dump(highly_variable_genes_data, f)

    # PCA
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,
        save=".png")
    pca_variance_ratio_data = {
        'pca_variance_ratio': adata.uns['pca']['variance_ratio'].tolist(),
    }
    with open(output_folder+'/pca_variance_ratio_data.json', 'w') as f:
        json.dump(pca_variance_ratio_data, f)

    sc.pl.pca(
        adata,
        color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
        save=".png"
    )

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(
        adata,
        color="sample",
        # Setting a smaller point size to get prevent overlap
        size=2,
        save=".png"
    )

    # Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
    sc.pl.umap(adata, color=["leiden"],
        save="_leiden.png")

    sc.pl.umap(
        adata,
        color=["leiden", "predicted_doublet", "doublet_score"],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
        save="_leiden2.png"
    )

    sc.pl.umap(
        adata,
        color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2,
        save="_leiden3.png"
    )

    # Extract UMAP coordinates
    umap_data = adata.obsm['X_umap']

    # Create a list of dictionaries for each cell, containing its UMAP coordinates and any other metadata you want to include
    cells = []

    if 'leiden' in adata.obs:  # Example of including cluster labels (assuming leiden clustering was done)
        clusters = adata.obs['leiden'].astype(str).tolist()
    else:
        clusters = ['None'] * len(umap_data)

    samples = adata.obs['sample'].astype(str).tolist()

    for i, (x, y) in enumerate(umap_data):
        cell_info = {
            'id': i,
            'x': float(x),
            'y': float(y),
            'cluster': clusters[i],  # Add any other metadata you need here
            'sample': samples[i]     # Extract the sample from which the data is taken
        }
        cells.append(cell_info)

    output = {
        "cells": cells,
        "type": "umap"
    }

    # Save the UMAP data as a JSON file
    with open(output_folder+'/umap_data.json', 'w') as f:
        json.dump(output, f)


    adata.write_h5ad(
        file.replace("data.h5", "scanpy/out.h5"),
        compression=hdf5plugin.FILTERS["zstd"],
        compression_opts=hdf5plugin.Zstd(clevel=5).filter_options
    )