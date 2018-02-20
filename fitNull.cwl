#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Run null model generation in genesis
hints:
  DockerRequirement:
    dockerPull: tmajarian/single-var@sha256:cfd43db702b6f067c21a54e3a4750ca341a3b11e2b81109c55b25752e6c98c4a
baseCommand: Rscript
inputs:
  script:
    type: string
    inputBinding:
      position: 1

  genotype_file:
    type: File
    inputBinding:
      position: 2
  
  phenotype_file:
    type: File
    inputBinding:
      position: 3

  outcome_name:
    type: string
    inputBinding:
      position: 4

  outcome_type:
    type: string
    inputBinding:
      position: 5

  covariates_string:
    type: string
    inputBinding:
      position: 6
      default: "NA"

  conditional_string:
    type: string
    inputBinding:
      position: 7
      default: "NA"

  ivars_string:
    type: string
    inputBinding:
      position: 8
      default: "NA"

  group_var:
    type: string
    inputBinding:
      position: 9
      default: "NA"

  sample_file:
    type: File
    inputBinding:
      position: 10
      default: "NA"

  label:
    type: string
    inputBinding:
      position: 11

  kinship_matrix:
    type: File
    inputBinding:
      prefix:
      position: 12

  id_col:
    type: string
    inputBinding:
      prefix:
      position: 13

outputs:
  model:
    type: File
    outputBinding:
      glob: "*.RDa"
