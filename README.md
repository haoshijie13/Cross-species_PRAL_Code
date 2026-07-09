# Cross-species PRAL Code

<img width="820" height="939" alt="image" src="https://github.com/user-attachments/assets/d5682a0d-62fe-4a68-a61d-2648ef443991" />

Code and analysis notebooks for cross-species transcriptomic and spatial transcriptomic study of vertebrate pallial organization, with a focus on reconstructing the primary sensory–allocortical (Pr–Al) molecular axis across evolution.

## Project summary

### Evolutionary Origins of Molecular Programs Underlying Brain Circuitry

The vertebrate pallium harbors independently evolved structures that, nevertheless, support strikingly similar sensory and cognitive circuit architectures. The mechanisms and evolutionary timing driving the emergence of these parallel pallial circuits remain unelucidated. Here, we integrated spatial transcriptomic and single-nucleus RNA-seq datasets from eight representative vertebrate species spanning approximately 500 million years of evolution to reconstruct the evolutionary assembly of the primary sensory-allocortical (Pr-Al) molecular axis, a conserved cortical hierarchical axis defined in our prior work. We uncovered that Pr- and Al-like neuronal identities are deeply conserved across sarcopterygians, encompassing all tetrapods and lobe-finned fishes. Intriguingly, this ancestral neuronal homology is uncoupled from spatially partitioned patterning: only tetrapods further compartmentalized these Pr- and Al- neurons into distinct pallium regions. Functional enrichment of Pr- and Al- gene programs uncovered a conserved tetrapod genetic core suite including MAPK signaling and axon guidance pathways. This core toolkit underwent sequential functional refinement from amphibians through mammals. Notably, mammals and birds convergently evolved association-cortical molecular profiles enriched for synaptic regulatory genes to support advanced cognitive functions. Collectively, this work delineates a stepwise vertebrate pallium evolutionary paradigm that explains how conserved molecular modules shape spatially organized brain circuitry across deep evolutionary time.

## Overview

This repository contains the main analysis notebooks and custom R helper functions used for:

- orthogroup-based gene harmonization across species
- integrated single-nucleus / single-cell RNA-seq analysis
- excitatory, inhibitory, and non-neuronal lineage-specific integration
- spatial transcriptomic region definition across multiple vertebrate species
- Pr–Al molecular axis analysis in RNA and spatial data
- RCTD-based spatial deconvolution / cell type transfer
- Mfuzz-based gene program dynamics analysis
- GO / KEGG functional annotation and enrichment analysis
- association-cortex related comparative analyses

The code is organized as a notebook-driven analysis workflow plus reusable plotting, smoothing, statistics, and enrichment utilities in `R_function/`.

## Repository structure

```text
MultiSpeciesCode/
├── Main_Code/
│   ├── 01.Orthgroup.ipynb
│   ├── 02.1.RNA_integrated.ipynb
│   ├── 03.1.RNA_integrated_EX.ipynb
│   ├── 03.2.RNA_integrated_IN.ipynb
│   ├── 03.3.RNA_integrated_NN.ipynb
│   ├── 04.1.ST_axol_lung_region_defined.ipynb
│   ├── 04.2.ST_bird_region_defined.ipynb
│   ├── 04.3.ST_fish_region_defined.ipynb
│   ├── 04.4.ST_lamp_region_defined.ipynb
│   ├── 04.5.ST_turt_region_defined.ipynb
│   ├── 04.6.ST_marm_mous_region_defined.ipynb
│   ├── 05.RNA_PrAl_analysis.ipynb
│   ├── 06.1.ST_PrAl_Gene.ipynb
│   ├── 06.2.ST_PrAl_RCTD.ipynb
│   ├── 06.3.ST_PrAl_RCTD_AllCelltype.ipynb
│   ├── 07.Gene_PrAl_analysis_mfuzz.ipynb
│   ├── 08.1.Prepared_GOKEGG.ipynb
│   ├── 08.2.GO_Anlysis.ipynb
│   └── 09.Association.ipynb
├── R_function/
│   ├── MultiSpecies.R
│   ├── Run_GO_KEGG.R
│   ├── ST_plot.R
│   ├── Seurat_function.R
│   ├── calculate_function.R
│   ├── flatmap_plot.R
│   ├── knn_function.R
│   └── plot_function.R
├── README_raw.md
└── README.md
```

## Analysis workflow

### 1. Orthogroup preparation
`01.Orthgroup.ipynb`

Builds orthogroup-based gene mapping tables from OrthoFinder output and prepares cross-species gene correspondence tables for downstream integration.

### 2. Multi-species RNA integration
`02.1.RNA_integrated.ipynb`

Performs integrated transcriptomic analysis across species using Seurat/Harmony-based workflows and builds the main integrated RNA reference.

### 3. Lineage-resolved RNA integration
- `03.1.RNA_integrated_EX.ipynb`: excitatory neuron integration
- `03.2.RNA_integrated_IN.ipynb`: inhibitory neuron integration
- `03.3.RNA_integrated_NN.ipynb`: non-neuronal cell integration

These notebooks further refine cross-species integration for major cellular compartments.

### 4. Spatial transcriptomic region definition
- `04.1.ST_axol_lung_region_defined.ipynb`
- `04.2.ST_bird_region_defined.ipynb`
- `04.3.ST_fish_region_defined.ipynb`
- `04.4.ST_lamp_region_defined.ipynb`
- `04.5.ST_turt_region_defined.ipynb`
- `04.6.ST_marm_mous_region_defined.ipynb`

These notebooks define spatial regions and region-level molecular organization in axolotl/lungfish, bird, fish, lamprey, turtle, and mammalian samples.

### 5. Pr–Al axis analysis in RNA data
`05.RNA_PrAl_analysis.ipynb`

Analyzes Pr–Al-related transcriptional structure in integrated RNA datasets.

### 6. Pr–Al axis analysis in spatial data
- `06.1.ST_PrAl_Gene.ipynb`: gene-level spatial Pr–Al analysis
- `06.2.ST_PrAl_RCTD.ipynb`: RCTD-based projection of Pr–Al identities into spatial data
- `06.3.ST_PrAl_RCTD_AllCelltype.ipynb`: extended RCTD/cell-type level analysis

### 7. Dynamic gene program analysis
`07.Gene_PrAl_analysis_mfuzz.ipynb`

Uses Mfuzz and related comparative analyses to study gene program dynamics along the Pr–Al axis.

### 8. Functional annotation and enrichment
- `08.1.Prepared_GOKEGG.ipynb`: prepares GO/KEGG annotation tables
- `08.2.GO_Anlysis.ipynb`: performs functional enrichment analysis and plotting

### 9. Association-cortex related comparison
`09.Association.ipynb`

Performs comparative analysis related to association-like molecular signatures and region/cell-type correspondence.

## R helper functions

The `R_function/` directory contains reusable helper functions used throughout the notebooks.

### `ST_plot.R`
Spatial plotting utilities, including:
- feature plotting in tissue coordinates
- metadata/column plotting
- figure export helpers
- mask-aware plotting support

### `knn_function.R`
kNN-based smoothing and projection utilities, including:
- iterative smoothing (`smooth_kNN`, `smooth_kNN2`)
- winner-take-all kNN assignment
- MAGIC-like smoothing
- ROC-style neighborhood evaluation
- object-level smoothing helpers

### `calculate_function.R`
General matrix/statistical utilities, including:
- group-wise matrix summarization
- scaling and normalization
- correlation helpers
- ordering and set conversion helpers
- z-test implementation

### `MultiSpecies.R`
Cross-species spatial plotting and enrichment helpers, including:
- multi-species spatial assay visualization
- multi-species metadata visualization
- custom enrichment plotting utilities

### `Run_GO_KEGG.R`
GO/KEGG analysis functions, including:
- GO enrichment
- KEGG enrichment
- gene ID conversion
- enrichment result plotting

### `flatmap_plot.R`
Flatmap visualization functions for marmoset, macaque, human, and mouse regional projections.

### `plot_function.R`
General visualization helpers, including quasirandom plots and correlation raster/dot plots.

### `Seurat_function.R`
Small utility functions for saving and restoring selected Seurat object slots/metadata.

## Main software dependencies

Based on the notebooks, the analysis primarily relies on the following R packages:

- Seurat
- harmony
- spacexr
- Mfuzz
- ggplot2
- patchwork
- ComplexHeatmap
- circlize
- dplyr
- Matrix
- ggsci
- scatterpie
- ggrastr
- scattermore
- RANN
- phangorn
- phytools
- ggtree
- AnnotationDbi
- GO.db
- org.Mm.eg.db
- org.Hs.eg.db
- UpSetR

Some notebooks also use species-specific annotations and external flatmap / shape resources.

## Data requirements and reproducibility notes

This repository is a code snapshot rather than a fully self-contained reproducible package.

Important points:

1. Many notebooks read data from absolute local paths such as `/mnt/gandan/...`, so the raw and intermediate data objects are not included here.
2. Several notebooks still call helper scripts that are not present in the current repository snapshot, including examples such as:
   - `segment_plot.R`
   - `Pr_Al.R`
   - `mutimodal.R`
3. Some notebook source references have been normalized during repository cleanup, but the missing helper scripts above are still required for full reruns.

Because of these external dependencies, the notebooks should be interpreted as analysis records and reusable code references unless the missing helper files and original datasets are restored.

## Suggested execution order

If you want to follow the main logic of the project, a reasonable reading/rerun order is:

1. `01.Orthgroup.ipynb`
2. `02.1.RNA_integrated.ipynb`
3. `03.*.RNA_integrated_*.ipynb`
4. `04.*.ST_*_region_defined.ipynb`
5. `05.RNA_PrAl_analysis.ipynb`
6. `06.*.ST_PrAl_*.ipynb`
7. `07.Gene_PrAl_analysis_mfuzz.ipynb`
8. `08.*.ipynb`
9. `09.Association.ipynb`

## Summary

This repository documents a multi-species comparative framework combining orthogroup harmonization, integrated RNA analysis, spatial mapping, Pr–Al axis reconstruction, and functional interpretation across vertebrate brains. The notebooks capture the full analytical storyline, while the `R_function/` folder provides the reusable building blocks for visualization, smoothing, enrichment, and comparative analysis.
