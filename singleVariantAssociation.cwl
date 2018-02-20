#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
inputs:
  this_genotype_files: File
  this_phenotype_file: File
  this_outcome_name: string
  this_outcome_type: string
  this_covariates_string: string
  this_conditional_string: string
  this_ivars_string: string
  this_group_var: string
  this_sample_file: File
  this_label: string
  this_kinship_matrix: File
  this_id_col: string
  this_test: string
  this_mac: int
  this_pval: string
  this_pval_threshold: float

outputs:
  pheno_with_conditional:
    type: File
    outputSource: conditionalPhenotype/new_phenotype_file
  alt_ref_conditional:
    type: File
    outputSource: conditionalPhenotype/alt_ref
  null_model:
    type: File
    outputSource: runFitNull/model
  all_associations:
    type: File
    outputSource: runSummary/allassoccsv
  top_associations: 
    type: File
    outputSource: runSummary/topassoccsv
  all_plots:
    type: File
    outputSource: runSummary/plots
  

steps:
  conditionalPhenotype:
    run: conditionalPhenotype.cwl
    in:
      genotype_files: this_genotype_files
      phenotype_file: this_phenotype_file
      sample_file: this_sample_file
      snps: this_conditional_string
      label: this_label
      script: "preprocess_conditional.R"
    out: [new_phenotype_file, alt_ref]


  runFitNull:
    run: fitNull.cwl
    in:
      genotype_file: this_genotype_files
      phenotype_file: conditionalPhenotype/new_phenotype_file
      outcome_name: this_outcome_name
      outcome_type: this_outcome_type
      covariates_string: this_covariates_string
      conditional_string: this_conditional_string
      ivars_string: this_ivars_string
      group_var: this_group_var
      sample_file: this_sample_file
      label: this_label
      kinship_matrix: this_kinship_matrix
      id_col: this_id_col
      script: "genesis_nullmodel.R"
    out: [model]

  runAssoc:
    run: assocTest.cwl
    in:
      gds_file: this_genotype_files
      null_file: runFitNull/model
      label: this_label
      test: this_test
      ivars_string: this_ivars_string
      mac: this_mac
      script: "association.R"
    out: [assoc_file]

  
  runSummary:
    run: summary.cwl
    in:
      pval: this_pval
      pval_threshold: this_pval_threshold
      label: this_label
      assoc: runAssoc/assoc_file
      script: "summary.R"
    out: [allassoccsv,topassoccsv,plots]