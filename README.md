# **Tertiary lymphoid structures harbor stem-like tumor-specific T cells**

**Alexander B. Afeyan†<sup>1**, Adi Nagler†, **Chloe R. Tu†<sup>2**, Gabriel Roberti De Oliveira, Berkay Simsek, Maxwell D. Seager, Nourhan El Ahmar, Haley E. Sax, Emma Lin, Amit Sud, Mehdi Borji, Cleo Forman, Sophia Liu, Patrick A. Ott, Toni K. Choueiri, Jennifer G. Abelin, Richard Burack, Shuqiang Li, Kenneth J. Livak, Svitlana Tyekucheva, Derin B. Keskin, Fei Chen, Michael B. Atkins, Jeremy M. Simon, Sabina Signoretti, Giacomo Oliveira°, David A. Braun‡°, Catherine J. Wu‡°

**Single-cell and bulk-RNAseq analysis<sup>1**, **spatial analysis<sup>2**

## Abstract
Tertiary lymphoid structures (TLS) are associated with improved response to immune checkpoint blockade across solid tumors1,2, but how they impact the phenotypic properties of tumor-specific T cells remains unknown. Across 24 treatment-naïve renal cell carcinoma (RCC) tumors, TLS-containing tumors were more heavily infiltrated by exhausted CD8+ T cells with a reduced terminal exhaustion transcriptional program compared to TLS- tumors. Specificity screening of 554 T cell clonotypes expanded within the microenvironment of 6 RCC tumors revealed 82 TCRs reactive against tumor cells and/or RCC antigens. A subset of tumor-specific T cell clonotypes (12%) was enriched within TLS, and these expressed an increased program of stem-like progenitor exhaustion, associated with favorable anti-tumor immunity. However, in 60 independent RCC tumors, macrophages within tumor margins of TLS-containing tumors bore an inferred immunosuppressive phenotype, and were colocalized with exhausted putative tumor-reactive T cells in a subgroup that was further analyzed, thus supporting this mode of immune evasion as a counterbalance to T cell immune pressure. Our data reveal TLS as reservoirs of tumor-specific T cells with stem-like progenitor features that could be leveraged by T cell immunotherapies.


## Contents
This repository contains all code used in the analyses of single-cell, bulk-RNA and spatial data associated with Afeyan Nagler Tu _et al_ (2025)

**Single-cell analysis as RMDs<sup>1**:
* Analysis of single-cell RNAseq data from Cohort 1 [available here](https://github.com/chloetutu/AfeyanNaglerTu_RCC/blob/main/single-cell/Afeyan%20Nagler%20Tu%20et%20al.%2C%202026%20single-cell%20processing%20Braun%202021.Rmd)
* Analysis of single-cell RNAseq data from Cohort 2 [available here](https://github.com/chloetutu/AfeyanNaglerTu_RCC/blob/main/single-cell/Afeyan%20Nagler%20Tu%20et%20al.%2C%202026%20single-cell%20processing%20Braun%202025.Rmd)

**Bulk-RNAseq analysis as RMD<sup>1**:
* Analysis of the bulk-RNAseq data [available here](https://github.com/chloetutu/AfeyanNaglerTu_RCC/blob/main/bulk-RNA/Afeyan%20Nagler%20Tu%20et%20al.%2C%202026%20bulk%20RNA%20processing.Rmd)

**Spatial analysis books rendered via Quarto<sup>2**:
* Integration of TCR and RNA data [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/00_TCR_RNA_integration_book/)
* RCTD and LIANA analysis [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/01_RCTD_LIANA_book/)
* Per puck annotation of spatial areas [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/02_area_annotation_book/)
* Clone and phenotype enrichment in spatial compartments [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/03a_spatial_enrichment_book/)
* QC of clone and phenotype enrichment in spatial compartments [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/03b_spatial_enrichment_QC_book/)
* Spatial bias of clones to specific TLS [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/03c_spatial_bias_book/)
* Ligand receptor colocalization analysis [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/04_ligand_receptor_colocalization_book/)
* Visualizations of spatial cell type features [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/05_plot_spatial_features_book/index.html)
* Visualizations of gene expression [available here](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/06_gene_expression_book/)
* Phenotype enrichment of clones in the TLS vs outside the TLS [available here
](https://chloetutu.github.io/AfeyanNaglerTu_RCC/spatial/07_shared_TCR_phenotypes_book)

Additionally, source quarto markdown documents (*.qmd) are supplied here under each book's subdirectory, e.g.:
https://github.com/chloetutu/AfeyanNaglerTu_RCC/tree/main/spatial/00_TCR_RNA_integration_book/source 
