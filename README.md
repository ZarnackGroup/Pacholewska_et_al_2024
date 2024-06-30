[![DOI](https://zenodo.org/badge/694054390.svg)](https://zenodo.org/doi/10.5281/zenodo.12597945)

# SF3B1 mutations primarily affect proximal alternative splice site in CLL and MDS patients.

## Abstract

Alternative splicing plays a critical role in generating transcriptome diversity, and aberrant splicing is frequently observed in cancer. Mutations in the splicing factor SF3B1 are particularly common in patients with chronic lymphocytic leukaemia (CLL) and myelodysplastic syndromes (MDS), with an opposite predictive model. In order to get insight into the effect of SF3B1 mutations on the transcriptome we used long-read sequencing data derived from MDS and CLL patient cohorts, as well as cell lines. Our results revealed that SF3B1 mutations specifically alter the usage of 3’ alternative splice sites especially within short proximity to the canonical splice sites. Moreover, disease-specific differences in the SF3B1 mutations effect seem to emerge from different transcriptomic profiles rather than different mechanism of SF3B1 mutation itself. Full isoform information enabled to predict the functional consequences of the aberrant splicing and gain mechanistic insights into the role of mutated SF3B1 in splicing.

## Source code

The source code is split into three parts:

* The analysis of the isoseq data ; Figures 1-4
* The analysis of the iCLIP data; Figure 5
* The integration of isoseq and iCLIP data; Figure 6

## Tools and packages

* The isoseq analysis for long read transcriptome sequencing was performed with the *IsoTools python module* ([GitHub](https://github.com/MatthiasLienhard/isotools))
* The iCLIP analysis was performed with the *BindingSiteFinder* R package ([GitHub](https://github.com/ZarnackGroup/BindingSiteFinder), [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/BindingSiteFinder.html)).

# Script structure

In this repository the Python notebooks, R-Scrips and Rmarkdown reports are provided. Scripts follow the outlined order:

## iCLIP data analysis

* 1 Binding site definition
  * Compute 5nt wide SF3B1 binding sites based on intial crosslink sites and iCLIP coverage
* 2 Binding site classification
  * Classify binding sites in of three bins according to their distance pattern
* 3 Coverage profile analysis
  * Compare SF3B1 and U2AF2 binding signals at selected genomic landmarks by direct coverage plots
* 4 Differential binding analysis
   * Analyse the binding site specific differences between the SF3B1 wt and mut condition

## Integration of isoseq and iCLIP

* 1 Splicing maps
  * Compute iCLIP binding maps across regulated 3'SS
* 2 Regulatory patterns
  * Analyse features of regulated and putative alternative 3'SS 


