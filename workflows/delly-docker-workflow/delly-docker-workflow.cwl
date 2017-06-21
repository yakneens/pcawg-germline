#!/usr/bin/env cwl-runner

class: CommandLineTool
id: bio_ppe-delly-Workflow
label: Finsen_delly-Workflow
cwlVersion: v1.0

doc: |
    ![build_status](https://img.shields.io/docker/build/weischenfeldt/bric-embl_delly_workflow.svg)
    ![docker_pulls](https://img.shields.io/docker/pulls/weischenfeldt/bric-embl_delly_workflow.svg)
    ![docker_builds](https://img.shields.io/docker/automated/weischenfeldt/bric-embl_delly_workflow.svg)

    A Dockstire version of the DELLY Structural Variant workflow used by the PPCG project.
    This is an updated and improved version of the PCAWG DELLY Workflow

    ### DELLY

    [DELLY](https://github.com/tobiasrausch/delly) is an integrated structural variant prediction method that can detect deletions, tandem duplications, inversions, insertions and translocations at single-nucleotide resolution in short-read massively parallel sequencing data.
    It uses paired-ends and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome.

dct:creator:
  '@id': http://orcid.org/0000-0003-3684-2659
  foaf:name: Francesco Favero
  foaf:mbox: francesco.favero@bric.ku.dk

dct:contributor:
  foaf:name: Etsehiwot Girum Girma
  foaf:mbox: Etsehiwot@finsenlab.dk

requirements:
  - class: DockerRequirement
    dockerPull: registry.hub.docker.com/weischenfeldt/finsen_delly_workflow:2.0.0_ppcg_test

hints:
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 16384
    outdirMin: 512000

inputs:
  run-id:
    type: string
    inputBinding:
      position: 1
      prefix: --run-id
  normal-bam:
    type: File
    inputBinding:
      position: 2
      prefix: --normal-bam
  tumor-bam:
    type: File
    inputBinding:
      position: 3
      prefix: --tumor-bam
  reference-gz:
    type: File
    inputBinding:
      position: 4
      prefix: --reference-gz
  reference-gc:
    type: File
    inputBinding:
      position: 5
      prefix: --reference-gc
  exclude-reg:
    type: File
    inputBinding:
      position: 6
      prefix: --exclude-reg
  sv-collection:
    type: File
    inputBinding:
      position: 7
      prefix: --sv-collection
  gc_wig:
    type: ["null", File]
    inputBinding:
      position: 8
      prefix: --gc_wig
  bin:
    type: ["null", int]
    inputBinding:
      position: 9
      prefix: --bin
      default: 200
  mem:
    type: ["null", int]
    inputBinding:
      position: 10
      prefix: --mem
  ncpu:
    type: ["null", int]
    inputBinding:
      position: 11
      prefix: --ncpu

outputs:
  bndout:
    type:
      type: array
      items: File
    outputBinding:
      glob: "bndout/*"
  highconfout:
    type:
      type: array
      items: File
    outputBinding:
      glob: "highconfout/*"
  sequenzaout:
    type:
      type: array
      items: File
    outputBinding:
      glob: "sequenzaout/*.tar.gz"
  dellyout:
    type:
      type: array
      items: File
    outputBinding:
      glob: "dellyout/*.tar.gz"
  archives:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.tar.gz"

baseCommand: [/usr/bin/start.sh, /usr/bin/launch_env.sh]