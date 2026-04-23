# Data and code for the "Human Nasal Microbiome Synthetic Communities" paper

The goal of this project was to study the assembly of Synthetic Communities (SynComs) that ressembled the Human Nasal Microbiome and study the metabolic interactions occuring between the bacteria in these SynComs.

## Info about the paper

**Abstract**:
<p align="justify">
The human nasal microbiome is a low-diversity ecosystem whose assembly principles and mechanisms of colonization resistance remain poorly understood. In particular, <em>Staphylococcus aureus</em> colonization shows high heterogeneity across individuals. We hypothesized that nutritional competition, strain-level diversity, and nutrient availability shape community stability and the ability of commensal species to inhibit <em>S. aureus</em>. To test this, we constructed 50 defined synthetic communities composed of representative human nasal bacteria differing in strain and species composition and tracked their temporal dynamics, metabolic profiles, and nutritional interactions. The composition of synthetic communities with 5-10 species showed robust and reproducible dynamics and converged in one of three stable states. Synthetic communities dominated by a specific strain of <em>Corynebacterium propinquum</em> were highly stable and consistently excluded <em>S. aureus</em>. Growth and coculture assays showed that <em>C. propinquum</em> outcompetes <em>S. aureus</em> under poor nutritional conditions resembling the nasal environment, whereas <em>S. aureus</em> dominates in nutrient-rich conditions. Metabolomics analyses revealed that nutritional competition, including siderophores utilization and amino acid limitation, likely underlies this colonization resistance. These results establish a tractable synthetic community model for the human nasal microbiome and identify nutrient-dependent competition and microbial metabolite production as key drivers of community structure and pathogen exclusion.
</p>

## Contact information

Marcelo Navarro-Diaz (marcelo.n.d@ciencias.unam.mx), Hannes Link (hannes.link@uni-tuebingen.de)

## Overview

<p align="justify">This repository contains the data and code for the analysis of the paper. There are two main folders: "Code" and "Data". The "Code" folder contains scripts running sequencing quality control and taxonomic assignment using Emu and subsequent statistical analyses. The "Data" folder contains all processed data.</p>

The raw data (.fastq files) is deposited in the Sequence Read Archive (SRA) of the NCBI under Bioproject accession number: PRJNA1370791.

The raw metabolomics data is deposited in MassIVE under accesion number: MSV000096776.

## Repository layout

**1. Code**

- **emu.sh**: Contains bash commands to run the sequences quality control and taxonomic assignation. 
- **code.R**: Has the R code necessary for the subsequent analysis, organized according to the figures of the paper. It is necessary to install the `nasalSynComsPkg` that contains helper functions for performing the analyses.
- **databases**: Contains the databases necessary to run emu. "16s_copies.csv" that contains the rRNA operon copy number for each species analysed in the sequencing data and "nose_sc_db_200824" that contains the emu database files for taxonomy assignment.

The recommended way to install `nasalSynComsPkg` is via **pak**, which automatically resolves
CRAN and Bioconductor dependencies:

```r
install.packages("pak")

pak::pak("github::marcelo-nd/nasalSynComsPkg")
```

### System requirements
- Windows: Rtools
- macOS: Xcode Command Line Tools

**2. Data**

- **Supplementary_Table_S1_HMP_ASV_table.biom**: ASV table containing diversity results for the Human Microbiome project data of Nasal Cavities (biom format).

- **Supplementary_Table_S1_HMP_ASV_table.csv**: ASV table containing diversity results for the Human Microbiome project data of Nasal Cavities (csv format).

- **Supplementary_Table_S2_Syncom_Inocula.xlsx**: Table containing the species/strains inoculated in each SynCom.

- **Supplementary Table S3. Growth Conditions.docx**: Table containing the growth conditions for strain activation from glycerol stocks.

- **Supplementary_Table_S4_Screening_OTU_table.csv**: OTU table containing diversity results for the screening of 50 SynComs.

- **Supplementary_Table_S5_Timepoints_OTU_table.csv**: OTU table containing diversity results for all time points for 20 SynComs.

- **Supplementary_Table_S6_SynCom_timepoints_metadata.csv**: Metadata table for time points samples for 20 SynComs.

- **Supplementary_Table_S7_Repetition_syncoms_OTU_table.csv**: Otu table containing diversity results for 5 SynComs.

- **Supplementary_Table_S8_Cocultures_OTU_table.csv**: OTU table containing divesity resutls for 3 cocultures in 3 SNM3, SNM10 and BHI.

- **Supplementary_Table_S9_aerobic_gcs.csv**: Table containing data for growth curves for aerobic strains.

- **Supplementary_Table_S10_anaerobic_gcs.csv**: Table containing data for growth curves for anaerobic strains.

- **Supplementary_Table_S11_Untargeted_feature_table.csv**: Untargeted metabolomics quant table for 20 SynComs.

- **Supplementary_Table_S12_Sirius_annotations.csv**: Sirius annotrations for untargeted metabolomics data for 20 SynComs.

- **Supplementary_Table_S13_Targeted_metabolomics_feature_table.xlsx**: Features table of metabolites for 5 SynComs.

- **Supplementary_Table_S14_C_propinquum_antiSMASH_BGCs.csv**: Table containing the predicted BGCs for the three <em>Corynebacterium propinquum</em> strains used in this study.




