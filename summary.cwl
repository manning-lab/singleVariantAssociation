#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Run summary after assoc test
hints:
  DockerRequirement:
    dockerPull: robbyjo/r-mkl-bioconductor:3.4.1
baseCommand: Rscript
inputs:
  script:
    type: string
    inputBinding:
      position: 1

  pval:
    type: string
    inputBinding:
      position: 2
  
  pval_threshold:
    type: float
    inputBinding:
      position: 3

  label:
    type: string
    inputBinding:
      position: 4

  assoc:
    type: File
    inputBinding:
      position: 5

outputs:
  allassoccsv:
    type: File
    outputBinding:
      glob: "*.assoc.csv"

  topassoccsv:
    type: File
    outputBinding:
      glob: "*.topassoc.csv"

  plots:
    type: File
    outputBinding:
      glob: "*.pdf"
