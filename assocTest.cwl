#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Run genesis association testing
hints:
  DockerRequirement:
    dockerPull: tmajarian/single-var@sha256:cfd43db702b6f067c21a54e3a4750ca341a3b11e2b81109c55b25752e6c98c4a
baseCommand: Rscript
inputs:
  script:
    type: string
    inputBinding:
      position: 1

  gds_file:
    type: File
    inputBinding:
      position: 2
  
  null_file:
    type: File
    inputBinding:
      position: 3

  label:
    type: string
    inputBinding:
      position: 4

  test:
    type: string
    inputBinding:
      position: 5
      default: "Score"

  ivars_string:
    type: string
    inputBinding:
      position: 6
      default: "NA"

  mac:
    type: int
    inputBinding:
      position: 7
      default: "5"

outputs:
  assoc_file:
    type: File
    outputBinding:
      glob: "*.assoc.RData"
