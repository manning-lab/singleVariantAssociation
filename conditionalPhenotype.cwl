#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Run summary after assoc test
hints:
  DockerRequirement:
    dockerPull: tmajarian/single-var@sha256:0d1369d43cdfe6ce762807cc475ab5dc530699d87c9c652b6744ad5f06f742a4
baseCommand: Rscript
inputs:

  script:
    type: string
    inputBinding:
      position: 1

  genotype_files:
    type: File[]
    inputBinding:
      position: 2
      itemSeparator: ","
      separate: false
  
  phenotype_file:
    type: File
    inputBinding:
      position: 3

  sample_file:
    type: string
    inputBinding:
      position: 4
      default: "NA"

  snps:
    type: string
    inputBinding:
      position: 5
      default: "NA"

  label:
  	type: string
  	inputBinding:
  		position: 6

outputs:
  new_phenotype_file:
    type: File
    outputBinding:
      glob: "*_phenotypes.csv"

  alt_ref:
    type: File
    outputBinding:
      glob: "*_alleles.txt"
