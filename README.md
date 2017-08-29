# SingleCellCrohns
scRNA seq of lymphocytes in Crohns disease

In this repository you will find the scripts for analysis of our single cell seq data using Seurat (v 1.4.0.13 Macosko 2015) and SCDE (v 1.99.4 Kharchenko 2016)
How to install as well as tutorials on Seurat can be found at: http://satijalab.org/seurat/
How to install as well as tutorials on SCDE can be found at: http://hms-dbmi.github.io/scde/

1.preparing_files_for_seurat is an R-script, with which a complete data matrix can be prepared
2.implement_seurat_and_filtering is an R-script, with which this data matrix is converted to Seurat-file and data is filtered for doublets and bad quality cells. Linear regression to regress out the effect of differences in number of UMI, percentage of Mitochondrial RNA and transcriptional differences per patient are corrected for.
3.JackStraw_analysis_cluster is an R-script in which JackStraw analysis is performed on a cluster computer in order to provide sufficient power
4.post_JackStraw is an R-script in which clusters are calculated based on JackStraw analyses
5.Change_ensemblgenenames_to_symbols is an R-script that describes how the ensembl gene names can be converted to gene symbols
6.Add_FACSdata is an R-script that describes how to add FACS data to the Seurat file
7.Cluster_Visualization is an R-script that describes various ways to visualize your results
8.DE_SCDE_blood_vs_mucosa is an R-script describing differential analysis of two subsets of data, this script can best be ran on a cluster computer, as it requires a lot of power
