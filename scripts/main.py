# https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
import scanpy as sc
import anndata as ad
import sys

if __name__=='__main__':
    folder = sys.argv[1]

    adatas = {}
    # Load the data
    sample_adata = sc.read_10x_mtx(folder, var_names='gene_symbols', cache=True)
    sample_adata.var_names_make_unique()
    adatas["sample"] = sample_adata
    
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()
    print(adata.obs["sample"].value_counts())
    print(adata)

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
    sc.settings.figdir = folder.replace("STAR/outSolo.out/Gene/filtered", "scanpy")
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=".png"
    )

    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save=".png")

    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.scrublet(adata, batch_key="sample")

    # Saving count data
    adata.layers["counts"] = adata.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    sc.pl.highly_variable_genes(adata,
        save=".png")

    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,
        save=".png")
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