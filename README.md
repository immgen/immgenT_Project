# immgenT Open Source Project

**A comprehensive molecular reference of T-cell states across the mouse immune system**

The **immgenT Open Source Project** (https://www.immgen.org/ImmGenT/) [is](https://www.immgen.org/ImmGenT/%29is) a collaborative effort of the Immunological Genome Project (ImmGen) to build a comprehensive reference atlas of mouse T cells across tissues, physiological states, infections, tumors, autoimmunity, and other immune perturbations.

The project profiled ~700,000 T cells from more than 700 samples using **single-cell RNA sequencing**, **128-plex CITE-seq surface-protein profiling**, and **paired TCRαβ sequencing**. Integration of these data defines a common molecular framework of major T-cell lineages and recurrent transcriptional states that can be compared across tissues and experimental conditions.

This repository contains code used for **processing, integration, analysis, visualization, and reference mapping of the immgenT dataset**, as well as analyses supporting the immgenT manuscripts.

## Repository organization

* `Data_processing_pipeline/` — processing and quality-control workflows for the immgenT single-cell data
* `Data_Integration/` — integration, dimensionality reduction, clustering, and reference construction
* `T-RBI/` — T-cell Reference-Based Integration tools for mapping external scRNA-seq datasets onto the immgenT reference
* `TCR_analysis/` — TCR repertoire analyses
* `GEX/` — gene-expression analyses
* `Cosmology_paper/` — analyses supporting the main immgenT reference manuscript
* `GP_paper/` — gene-program analyses
* `immgenT-CD4/`, `immgenT-CD8/`, `Treg_paper/` — lineage-specific analyses and companion manuscripts

Individual directories contain the scripts and workflows associated with each analysis.

## Data

Raw and processed immgenT sequencing data are available through the **Gene Expression Omnibus (GEO)**:

**GSE297097**

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297097

## Explore immgenT

The main immgenT portal provides access to the reference atlas and associated interactive resources:

**https://www.immgen.org/ImmGenT/**

Available tools include:

* **Rosetta2** — interactive exploration of RNA and CITE-seq data
* **T-RBI** — mapping and annotation of external scRNA-seq datasets in the immgenT reference
* **Skyline** — pseudobulk gene-expression profiles across immgenT cell states
* **TCR Browser** — exploration of TCR sequences across the immgenT dataset

## T-RBI: map your data to immgenT

T-RBI (**T-cell Reference-Based Integration**) enables users to project external single-cell RNA-seq datasets onto the fixed immgenT reference and obtain:

* immgenT lineage and cluster annotations
* annotation confidence scores
* discovery scores for cells poorly represented in the reference
* coordinates in the all-T and lineage-specific immgenT MDE maps

The public T-RBI interface is available at:

**https://rstats.immgen.org/immgenT_integration_analysis/**

## Citation

The immgenT reference and associated analyses are described in the immgenT manuscript series. Citation information will be updated here as the manuscripts are published.

For the original description of the project:

immgenT Cosmology: Magill et al. **immgenT: A Comprehensive Reference of Convergent T-cell States in the Mouse.** https://www.biorxiv.org/content/10.64898/2026.01.30.702892

immgenT TCR: Croze et al. **The αβTCR repertoire at scale in the immgenT dataset.** https://www.biorxiv.org/content/10.64898/2026.01.30.702900

immgenT CD8: Galletti et al. **The CD8 immgenT framework as a universal reference of mouse CD8αβ T cell differentiation states.** https://www.biorxiv.org/content/10.64898/2026.02.02.703365

immgenT Treg: Freuchet et al. **A Reference Landscape of Regulatory T Cell States in Mice.** https://www.biorxiv.org/content/10.64898/2026.01.30.702856

immgenT CD4: Mehrotra et al. **A Common Reference Landscape of Mouse CD4+ T Cells.** https://www.biorxiv.org/content/10.64898/2026.05.29.728553
