#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "cgpmap"

label: "CGP BWA-mem mapping flow"

cwlVersion: v1.0

doc: |
    ![build_status](https://quay.io/repository/wtsicgp/dockstore-cgpmap/status)
    A Docker container for the CGP BWA-mem mapping flow. See the [dockstore-cgpmap](https://github.com/cancerit/dockstore-cgpmap) website for more information.

dct:creator:
  "@id": "http://orcid.org/0000-0002-5634-1539"
  foaf:name: Keiran M Raine
  foaf:mbox: "keiranmraine@gmail.com"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/dockstore-cgpmap:1.0.6"

hints:
  - class: ResourceRequirement
    coresMin: 1 # works but long, 8 recommended
    ramMin: 15000 # good for WGS human ~30-60x
    outdirMin: 5000000 # unlikely any BAM processing would be possible in less

inputs:
  reference:
    type: File
    doc: "The core reference (fa, fai, dict) as tar.gz"
    inputBinding:
      prefix: -reference
      position: 1
      separate: true

  bwa_idx:
    type: File
    doc: "The BWA indexes in tar.gz"
    inputBinding:
      prefix: -bwa_idx
      position: 2
      separate: true

  sample:
    type: string
    doc: "Sample name to be included in output [B|CR]AM header, also used to name final file"
    inputBinding:
      prefix: -sample
      position: 3
      separate: true

  scramble:
    type: string?
    doc: "Options to pass to scramble when generating CRAM output, see scramble docs"
    default: ''
    inputBinding:
      prefix: -scramble
      position: 4
      separate: true
      shellQuote: true

  bwa:
    type: string?
    default: ' -Y -K 100000000'
    doc: "Mapping and output parameters to pass to BWA-mem, see BWA docs, default ' -Y -K 100000000'"
    inputBinding:
      prefix: -bwa
      position: 5
      separate: true
      shellQuote: false

  cram:
    type: boolean
    doc: "Set if output should be in CRAM format instead of BAM, see 'scramble' for tuning parameters."
    inputBinding:
      prefix: -cram
      position: 6

  bams_in:
    type:
    - 'null'
    - type: array
      items: File
    doc: "Can be BAM, CRAM, fastq (paired or interleaved), BAM/CRAM can be mixed together but not FASTQ."
    inputBinding:
      position: 7

outputs:
  out_bam:
    type: File
    outputBinding:
      glob: $(inputs.sample).bam

  out_bai:
    type: File
    outputBinding:
      glob: $(inputs.sample).bam.bai

  out_bas:
    type: File
    outputBinding:
      glob: $(inputs.sample).bam.bas

  out_md5:
    type: File
    outputBinding:
      glob: $(inputs.sample).bam.md5

  out_met:
    type: File
    outputBinding:
      glob: $(inputs.sample).bam.met

  out_maptime:
    type: File
    outputBinding:
      glob: $(inputs.sample).bam.maptime

baseCommand: ["/opt/wtsi-cgp/bin/ds-wrapper.pl"]
