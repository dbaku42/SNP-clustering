# nf-core/snp_clustering

**nf-core/snp_clustering** is a pipeline built using Nextflow that follows [nf-core](https://nf-co.re) best practices.  
It performs unsupervised clustering of individuals based on SNP genotypes using PCA, t-SNE, and K-Means.

---

## ✨ Key Features

- Accepts both compressed and uncompressed VCF files as input.
- Encodes genotypes numerically: `0`, `1`, `2`; `-1` for missing data.
- Filters variants based on quality (missing rate > 5% or low variance).
- Reduces dimensionality using PCA and standard scaling.
- Performs clustering with K-Means (includes elbow plot).
- Visualizes results with PCA and t-SNE scatterplots.

---

## 🔄 Workflow Summary

graph TD\
  A[VCF Input File] --> B[Genotype Matrix Extraction & Cleaning]\
  B --> C[PCA on Scaled Data]\
  C --> D[K-Means Clustering]\
  D --> F[PCA Cluster Plot]\
  D --> E[t-SNE Visualization]

## 📥 Input
A .vcf or .vcf.gz file containing diploid, bi-allelic SNPs.

Genotype encoding:

- 0: Homozygous reference

- 1: Heterozygous

- 2: Homozygous alternate

- -1: Missing data

### 📤 Output Files

| **File**                    | **Description**                                 |
|-----------------------------|--------------------------------------------------|
| `genotype_matrix.npy`       | Cleaned genotype matrix                         |
| `samples.txt`               | List of sample IDs                              |
| `variants_filtered.txt`     | List of filtered variant IDs                    |
| `pca_coords.csv`            | PCA coordinates of samples                      |
| `elbow_plot.png`            | Elbow plot for optimal `k` selection            |
| `cluster_assignments.csv`   | K-Means cluster labels                          |
| `pca_clustered_plot.png`    | PCA plot colored by cluster                     |
| `tsne_coords.csv`           | 2D coordinates from t-SNE                       |
| `tsne_clustered_plot.png`   | t-SNE plot colored by cluster                   |

### 🧱 Implementation Details

The pipeline is built using **Nextflow DSL2** and structured in modular `process` blocks. It relies on Python for all computation and visualization tasks.

**Python tools used:**

- `scikit-allel`: Genotype extraction from VCF files  
- `numpy`, `pandas`: Matrix transformation and data handling  
- `scikit-learn`: PCA, KMeans clustering, and t-SNE embedding  
- `matplotlib`: Plotting (PCA, t-SNE, elbow method)
### 🔧 Requirements

- **Nextflow** ≥ 21.04.0
- **Python** ≥ 3.8 with the following packages:
  - `numpy`
  - `pandas`
  - `scikit-allel`
  - `scikit-learn`
  - `matplotlib`

> ⚙️ Environments are managed via **Docker**, **Conda**, or **Singularity** profiles to ensure reproducibility and ease of setup.
## Author
**Donald Baku** Github: https://github.com/dbaku42
