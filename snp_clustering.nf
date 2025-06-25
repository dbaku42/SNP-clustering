#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define input parameters with default values
params.vcf_file = "input.vcf.gz" // Can be .vcf or .vcf.gz
params.out_dir = "results"
params.n_clusters = 4            // Desired number of clusters for K-Means
params.max_k_elbow = 10          // Max k to test for the elbow plot

/*
========================================================================================
    Process 1: Extract, filter, and prepare the genotype matrix
========================================================================================
*/
process extract_genotype_matrix {
    publishDir "${params.out_dir}/genotype_data", mode: 'copy'

    input:
    path vcf_file

    output:
    path "genotype_matrix.npy"   emit: genotype_matrix
    path "samples.txt"           emit: samples_list
    path "variants_filtered.txt" emit: variants_list

    script:
    """
    #!/usr/bin/env python
    import allel
    import numpy as np

    print(f"Processing VCF file: {vcf_file}")
    vcf = allel.read_vcf(f"{vcf_file}", fields=['calldata/GT', 'variants/ID', 'samples'])
    gt = allel.GenotypeArray(vcf['calldata/GT'])
    gm = gt.to_n_alt(1)
    print(f"Initial genotype matrix shape: {gm.shape}")

    # Filter 1: Drop SNPs with NA >= 5%
    missing_rate = np.count_nonzero(gm == -1, axis=1) / gm.shape[1]
    na_filter_mask = missing_rate < 0.05
    gm_filtered_na = gm[na_filter_mask]
    variants_filtered_na = vcf['variants/ID'][na_filter_mask]
    print(f"Shape after NA filter: {gm_filtered_na.shape}")

    # Filter 2: Drop SNPs with variance < 0.05
    gm_imputed_for_var = np.copy(gm_filtered_na).astype(float)
    for i in range(gm_imputed_for_var.shape[0]):
        variant_row = gm_imputed_for_var[i, :]
        non_missing = variant_row[variant_row != -1]
        if non_missing.size > 0:
            row_mean = non_missing.mean()
            variant_row[variant_row == -1] = row_mean
            gm_imputed_for_var[i, :] = variant_row
    
    variance = np.var(gm_imputed_for_var, axis=1)
    var_filter_mask = variance >= 0.05
    gm_final = gm_filtered_na[var_filter_mask]
    variants_final = variants_filtered_na[var_filter_mask]
    print(f"Shape after variance filter: {gm_final.shape}")

    np.save("genotype_matrix.npy", gm_final)

    samples = vcf['samples']
    with open("samples.txt", "w") as f:
        for sample_id in samples: f.write(f"{sample_id}\\n")

    with open("variants_filtered.txt", "w") as f:
        for variant_id in variants_final: f.write(f"{variant_id}\\n")
    """
}

/*
========================================================================================
    Process 2: Standard Scaling and Principal Component Analysis (PCA)
========================================================================================
*/
process perform_pca {
    publishDir "${params.out_dir}/pca_results", mode: 'copy'

    input:
    path genotype_matrix
    path samples_list

    output:
    path "pca_coords.csv" emit: pca_results_csv

    script:
    """
    #!/usr/bin/env python
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    gm = np.load(f"{genotype_matrix}")
    samples = pd.read_csv(f"{samples_list}", header=None)[0].tolist()
    
    gm_transposed = gm.T
    scaler = StandardScaler()
    gm_scaled = scaler.fit_transform(gm_transposed)
    
    # We run PCA for enough components to be used in clustering and plotting
    n_components_pca = min(10, len(samples) - 1, gm_scaled.shape[1])
    pca = PCA(n_components=n_components_pca)
    coords = pca.fit_transform(gm_scaled)
    
    explained_variance = pca.explained_variance_ratio_
    
    # Save results to a CSV file
    pc_cols = [f"PC{i+1}" for i in range(coords.shape[1])]
    coords_df = pd.DataFrame(coords, columns=pc_cols, index=samples)

    # Also save explained variance
    coords_df.to_csv("pca_coords.csv", index_label="SampleID")
    
    print("PCA calculation complete and coordinates saved.")
    """
}

/*
========================================================================================
    Process 3: K-Means Cluster Analysis
========================================================================================
*/
process perform_clustering {
    publishDir "${params.out_dir}/cluster_analysis", mode: 'copy'

    input:
    path pca_results_csv

    output:
    path "cluster_assignments.csv" emit: cluster_assignments_csv
    path "elbow_plot.png"          emit: elbow_plot
    path "pca_clustered_plot.png"  emit: pca_clustered_plot

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.cluster import KMeans

    pca_df = pd.read_csv(f"{pca_results_csv}", index_col="SampleID")

    # --- Elbow Method to find optimal K ---
    # Use the first 10 PCs or fewer if not available
    n_pcs_to_use = min(10, pca_df.shape[1])
    pca_data_for_clustering = pca_df.iloc[:, :n_pcs_to_use]
    
    inertia = []
    max_k = min(${params.max_k_elbow}, len(pca_df) - 1)
    k_range = range(1, max_k + 1)
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto').fit(pca_data_for_clustering)
        inertia.append(kmeans.inertia_)

    plt.figure(figsize=(10, 6))
    plt.plot(k_range, inertia, 'bo-')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Inertia (Within-cluster sum of squares)')
    plt.title('Elbow Method For Optimal k')
    plt.grid(True)
    plt.savefig("elbow_plot.png")
    plt.close()

    # --- Perform Final Clustering with specified n_clusters ---
    n_clusters = ${params.n_clusters}
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init='auto').fit(pca_data_for_clustering)
    clusters = kmeans.labels_

    # Save cluster assignments
    assignments_df = pd.DataFrame({'Cluster': clusters}, index=pca_df.index)
    assignments_df.to_csv("cluster_assignments.csv")

    # --- Create Clustered PCA Plot ---
    plt.figure(figsize=(12, 10))
    scatter = plt.scatter(pca_df['PC1'], pca_df['PC2'], c=clusters, cmap='viridis', alpha=0.8)
    plt.xlabel(f"PC1")
    plt.ylabel(f"PC2")
    plt.title(f'PCA colored by K-Means Clusters (k={n_clusters})')
    plt.legend(handles=scatter.legend_elements()[0], labels=[f'Cluster {i}' for i in range(n_clusters)], title="Clusters")
    plt.grid(True)
    plt.savefig("pca_clustered_plot.png")
    plt.close()
    """
}

/*
========================================================================================
    Process 4: t-SNE Visualization with Cluster Colors
========================================================================================
*/
process perform_tsne_visualization {
    publishDir "${params.out_dir}/tsne_results", mode: 'copy'

    input:
    path pca_results_csv
    path cluster_assignments_csv

    output:
    path "tsne_coords.csv"      emit: tsne_results_csv
    path "tsne_clustered_plot.png" emit: tsne_plot

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.manifold import TSNE

    pca_df = pd.read_csv(f"{pca_results_csv}", index_col="SampleID")
    cluster_df = pd.read_csv(f"{cluster_assignments_csv}", index_col="SampleID")
    
    # Align data before t-SNE
    data = pca_df.join(cluster_df)

    n_samples = len(data)
    perplexity_value = min(30.0, n_samples - 1.0)
    
    tsne = TSNE(n_components=2, perplexity=perplexity_value, n_iter=1000, learning_rate='auto', init='pca', random_state=42)
    tsne_coords = tsne.fit_transform(data.iloc[:, :-1]) # Use PCA coords for t-SNE

    tsne_df = pd.DataFrame(tsne_coords, columns=["tSNE1", "tSNE2"], index=data.index)
    tsne_df.to_csv("tsne_coords.csv")

    # Create plot colored by cluster assignment
    plt.figure(figsize=(12, 10))
    clusters = data['Cluster']
    n_clusters = len(clusters.unique())
    scatter = plt.scatter(tsne_df['tSNE1'], tsne_df['tSNE2'], c=clusters, cmap='viridis', alpha=0.8)
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.title('t-SNE colored by K-Means Clusters')
    plt.legend(handles=scatter.legend_elements()[0], labels=[f'Cluster {i}' for i in range(n_clusters)], title="Clusters")
    plt.grid(True)
    plt.savefig("tsne_clustered_plot.png")
    plt.close()
    """
}


/*
========================================================================================
    Workflow Definition
========================================================================================
*/
workflow {
    vcf_ch = file(params.vcf_file)

    // 1. Extract and filter genotype data
    extract_genotype_matrix(vcf_ch)

    // 2. Perform scaling and PCA
    perform_pca(
        extract_genotype_matrix.out.genotype_matrix,
        extract_genotype_matrix.out.samples_list
    )

    // 3. Perform K-Means clustering on PCA results
    perform_clustering(perform_pca.out.pca_results_csv)

    // 4. Perform t-SNE and visualize with cluster labels
    perform_tsne_visualization(
        perform_pca.out.pca_results_csv,
        perform_clustering.out.cluster_assignments_csv
    )
}
