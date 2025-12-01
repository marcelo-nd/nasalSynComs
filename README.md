# Data and code for the "Human Nasal Microbiome Synthetic Communities" project

The goal of this project was to study the assembly Synthetic Communities (SynComs) that ressembled the Human Nasal Microbiome and study the metabolic interactions occuring between the bacteria in these SynComs. 

## Info about the paper

**Abstract**:
The human nasal microbiome is a low-diversity ecosystem whose assembly principles and mechanisms of colonization resistance remain poorly understood. In particular, Staphylococcus aureus colonization shows high heterogeneity across individuals. We hypothesized that nutritional competition, strain-level diversity, and nutrient availability shape community stability and the ability of commensal species to inhibit S. aureus. To test this, we constructed 50 defined synthetic communities composed of representative human nasal bacteria differing in strain and species composition and tracked their temporal dynamics, metabolic profiles, and nutritional interactions. Community assembly followed reproducible trajectories characterized by three stable states in which the specific strain of Corynebacterium propinquum strongly determined whether this species dominated. Synthetic communities dominated by C. propinquum were highly stable and consistently excluded S. aureus. Growth and coculture assays showed that C. propinquum outcompetes S. aureus under nutrient-limited conditions resembling the nasal environment, whereas S. aureus prevails only in nutrient-rich conditions. Metabolomics analyses revealed that nutritional competition, including siderophores utilization and amino acid limitation, likely underlies this colonization resistance. These results establish a tractable synthetic community model for the human nasal microbiome and identify nutrient-dependent competition and microbial metabolite production as key drivers of community structure and pathogen exclusion.

## Contact information

Marcelo Navarro-Diaz (marcelo.n.d@ciencias.unam.mx)
Hannes Link (hannes.link@uni-tuebingen.de)

## Overview

This repository contains the data and code for the analysis of the paper. There are two main folders: Code and data. The first folder contains scripts for R and running sequencing quality control and taxonomic assignment using Emu. The second folder contains all processed data.

The raw data (.fastq files) is deposited in the SRA of the NCBI with Bioporject accession number: PRJNA1370791. 

## Repository layout

**1. Code**

- **emu.sh**: Contains bash commands to run the sequences quality control and taxonomic assignation. 
- **code.R**: Has the R code neccesary for the subsequent analysis, organized according to the figures of the paper. 

**2.Data**

- **1_screening_otu_table.csv**: Otu table containing diversity results for the screening of 50 SynComs.

- **2_nasal_syncom_strains.xlsx**: Table containing the species/strains contained in each SynCom.

- **3_timepoints_otu_table.csv**: Otu table containing diversity results for all time points for 20 SynComs.

- **4_timepoints_metadata.csv**: Metadata table for time points samples for 20 SynComs.

- **5_untargeted_quant_table.csv**: Untargeted metabolomics quant table for 20 SynComs.

- **6_sirius_annotations.csv**: Sirius annotrations for untargeted metabolomics data for 20 SynComs.

- **7_repetition_syncoms_otu_table.csv**: Otu table containing diversity results for 5 SynComs.

- **8_20251030_12C_Nasal_targeted_metabolomics_data_002.xlsx**: Features table of metabolites for 5 SynComs.

- **9_hmp_asv_table.biom**: ASV table containing diversity results for the Human Microbiome project data of Nasal Cavities.

- **10_cocultures_otu_table.csv**: OTU table containing divesity resutls for 3 cocultures in 3 SNM3, SNM10 and BHI.



