# https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
import scanpy as sc
import anndata as ad
import sys
import hdf5plugin
import json

if __name__=='__main__':
    folder = sys.argv[1]

    groups_file_path = folder.replace("STAR/outSolo.out/Gene/filtered", "groups.json")
    
    groups_data = []
    with open(groups_file_path) as f:
        groups_data = json.load(f)

    adatas = {}
    for group in groups_data:
        # print(group["name"])
        sample_folder = folder.replace("outSolo.out", group["name"]+"Solo.out")
        sample_adata = sc.read_10x_mtx(sample_folder, var_names='gene_symbols', cache=True)
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
    output_folder = folder.replace("STAR/outSolo.out/Gene/filtered", "scanpy")
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
        'pca_variance': adata.uns['pca']['variance'].tolist(),
        'pca_variance_ratio': adata.uns['pca']['variance_ratio'].tolist(),
        'n_pcs': 50
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
    pca_plot_data = {
        'sample': adata.obs['sample'].tolist(),
        'pct_counts_mt': adata.obs['pct_counts_mt'].tolist(),
        'pca_0': adata.obsm['X_pca'][:, 0].tolist(),
        'pca_1': adata.obsm['X_pca'][:, 1].tolist(),
        'pca_2': adata.obsm['X_pca'][:, 2].tolist(),
        'pca_3': adata.obsm['X_pca'][:, 3].tolist()
    }
    with open(output_folder+'/pca_plot_data.json', 'w') as f:
        json.dump(pca_plot_data, f)

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

    adata.write_h5ad(
        folder.replace("STAR/outSolo.out/Gene/filtered", "scanpy/out.h5"),
        compression=hdf5plugin.FILTERS["zstd"],
        compression_opts=hdf5plugin.Zstd(clevel=5).filter_options
    )