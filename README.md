# Single Variant Association -- Linear Mixed Models 

## Description 

This workflow performs a single variant association analysis of genotype data with a single phenotype using linear mixed models with fixed effects. The primary code is written in R using the GENESIS package for model fitting and association testing. The workflow can either generate a null model from phenotype and relatedness data or use a pregenerated null model.

### Authors

This workflow is produced and maintained by the [Manning Lab](https://manning-lab.github.io/). Contributing authors include:

* Tim Majarian (tmajaria@broadinstitute.org)
* Jasen Jackson (jasenjackson97@gmail.com)
* Alisa Manning (amanning@broadinstitute.org).

## Dependencies

### Workflow execution

* [WDL](https://software.broadinstitute.org/wdl/documentation/quickstart)
* [Chromwell](http://cromwell.readthedocs.io/en/develop/)

### R packages

* [GENESIS](https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html)
* [GWASTools](https://www.bioconductor.org/packages/release/bioc/html/GWASTools.html)
* [SeqArray](https://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)
* [SeqVarTools](https://www.bioconductor.org/packages/release/bioc/html/SeqVarTools.html)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [qqman](https://cran.r-project.org/web/packages/qqman/index.html)

## Main Functions

### preprocessConditional
This function preprocesses the phenotype file for a conditional analysis on some snps, extracting the genotype information and adding to the phenotypes for each sample.

Inputs:
* genotype_files : comma separated list of gds file paths, this list must contain gds files for every snp
* phenotype_file : phenotype file to be edited (file, CSV/TSV)
* sample_file : text file with list of sample ids to include, one per line
* snps : comma separated list with the form chr_number:position,chr_number:position (string)

Outputs :
* new_phenotype_file : phenotype.data input with appended dosage columns for snps of interest (.csv)
* alt_ref : text file with alternate and reference alleles for each snp of interest (.txt)

### fitNull

This function generates a null model to be used in association testing in Genesis

Inputs:
* genotype_files : genotype data for all samples (array of VCF or GDS file)
* phenotype_file : phenotype data for all samples to be included in analysis (CSV or TSV file)
* outcome : the outcome to be tested (string)
* outcome_type : the type of outcome being tested (dichotomous or continuous)
* covariates_string : covariates to condition on in linear mixed modeling (comma separated string, default = '')
* sample_file : a file containing a list of sample ids (matching the genotype and phenotype files) to be included, one per line (.txt)
* label : prefix for output filename (string)
* kinship_file : relatedness measures for all samples (CSV or TSV file)
* id_col : column name of id column in phenotype file (string)

Outputs:
* model : generated null model (.RDa)

### assocTest

This function performs an association test to generate p-values for each variant included.

Inputs:
* gds_file : a genotype file containing data for all samples are variants to be tested (.gds)
* null_file : output of the function *fitNull* or a pregenerated null model (.RDa)
* label : prefix for output filename (string)
* test : statistical test (Score or Wald, default = Score)
* mac : minimum minor allele count for variants to be included in analysis (int, default = 5)

Outputs:
* assoc : an RData file of associations results (.RData)

### summary

Generate a summary of association results including quantile-quantile and manhattan plots for variants subseted by minor allele frequency (all variants, maf < 5%, maf >= 5%). Also generates CSV files of all and the top associated variants.

Inputs:
* pval : the p-value column in the output of assocTest, this should be the statistical test with ".pval" appended (string, Score -> Score.pval, Wald -> Wald.pval)
* pval_threshold : p-value threshold for the returning top associations, top association output will include only variants with a p-value less than the threshold (float, default = 0.0001)
* label : prefix for output filename (string)
* assoc : output of assocTest (Array[.RData])

## Other workflow inputs

* this_memory : amount of memory in GB for each execution of a task (int)
* this_disk : amount of disk space in GB to allot for each execution of a task (int)



