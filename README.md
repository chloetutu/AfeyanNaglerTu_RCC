# **Tertiary lymphoid structures are reservoirs of stem-like tumor-specific T cells in renal cell carcinoma**

**Alexander B. Afeyan†**, Adi Nagler†, **Chloe R. Tu†**, Gabriel Roberti De Oliveira, Berkay Simsek, Maxwell D. Seager, Nourhan El Ahmar, Haley E. Sax, Emma Lin, Amit Sud, Mehdi Borji, Cleo Forman, Sophia Liu, Patrick A. Ott, Toni K. Choueiri, Jennifer G. Abelin, Richard Burack, Shuqiang Li, Kenneth J. Livak, Svitlana Tyekucheva, Derin B. Keskin, Fei Chen, Michael B. Atkins, Jeremy M. Simon, Sabina Signoretti, Giacomo Oliveira°, David A. Braun‡°, Catherine J. Wu‡°

## Abstract
Tertiary lymphoid structures (TLS) are associated with improved response to immune checkpoint blockade across solid tumors1–3, but their impact on tumor-specific T cell responses remains incompletely defined. By tracking the spatial phenotype of T cell clones with confirmed anti-tumor specificity, we report the architecture of T cell responses in renal cell carcinoma (RCC). Across 24 RCC tumors, TLS-containing (TLS+) tumors were more infiltrated by exhausted CD8 T cells, but with a reduced terminal exhaustion transcriptional program compared to TLS- tumors. To investigate the effect of TLS localization on RCC-specific T cells, we screened 581 CD8+ and CD4+ T cell receptor (TCR) α/β pairs across six patients and identified 72 tumor-specific (TS) and 21 tumor-antigen-specific TCRs. Among these bona fide TS-TCR clonotypes, 45.8% of clonotypes infiltrated tumors, and 10.2% localized within TLS by spatial transcriptomics. TS-TCRs in TLS+ tumors preferentially colocalized with plasma cell transcripts (P<0.0001) compared to viral-reactive TCRs, and were associated with upregulated progenitor exhausted. However, in counterbalance, exhausted TCR clonotypes also colocalized with macrophages of inferred immunosuppressive phenotype concentrated within a myeloid boundary of the tumor (P=3x10-14). Our data reveal TLS as reservoirs of tumor-specific T cells with stem-like features that could be leveraged by T cell immunotherapies.

## Contents
This repository contains all code used in the analyses of single-cell, bulk-RNA and spatial data associated with Afeyan Nagler Tu _et al_ (2025)

Single-cell analysis as RMDs:
* Analysis of single-cell RNAseq data from Cohort 1 [available here](https://github.com/chloetutu/AfeyanNaglerTu_RCC/blob/main/single-cell/Afeyan%20Nagler%20Tu%20et%20al.%2C%202026%20single-cell%20processing%20Braun%202021.Rmd)
* Analysis of single-cell RNAseq data from Cohort 2 [available here](https://github.com/chloetutu/AfeyanNaglerTu_RCC/blob/main/single-cell/Afeyan%20Nagler%20Tu%20et%20al.%2C%202026%20single-cell%20processing%20Braun%202025.Rmd)

Bulk-RNAseq analysis as RMD:
* Analysis of the bulk-RNAseq data [available here](https://github.com/chloetutu/AfeyanNaglerTu_RCC/blob/main/bulk-RNA/Afeyan%20Nagler%20Tu%20et%20al.%2C%202026%20bulk%20RNA%20processing.Rmd)

Spatial analysis books rendered via Quarto:
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
