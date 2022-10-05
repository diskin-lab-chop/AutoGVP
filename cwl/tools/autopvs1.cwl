cwlVersion: v1.2
class: CommandLineTool
id: kfdrc-autopvs1
label: AutoPVS1 for VEP Output
doc: |
  "This is a limited application of running autopvs1 on an already annotated VEP file"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'migbro/autopvs1:v0.2.0a'

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      tar --use-compress-program="pigz -p 8" -xf
  - position: 1
    shellQuote: false
    valueFrom: >-
      && python3 /autopvs1/config_create.py
  - position: 2
    shellQuote: false
    valueFrom: >-
      && python3 /autopvs1/autoPVS1_from_VEP_vcf.py
  - position: 3
    shellQuote: false
    valueFrom: >-
      > $(inputs.output_basename).autopvs1.tsv

inputs:
  autopvs1_db: { type: File, doc: "Annovar Database with at minimum required resources to InterVar", inputBinding: { position: 0 }}
  autopvs1_db_str: { type: 'string?', doc: "Name of dir created when annovar db is un-tarred", inputBinding: { position: 1, prefix: "--data_dir" }, default: "data"}
  vep_vcf: { type: File, doc: "VEP annotated input file", inputBinding: { position: 2, prefix: "--vep_vcf"} }
  genome_version:  { type: ['null', { type: enum, symbols: ["hg38","GRCh38"], name: "genome_version" } ], doc: "Genome reference build version",
    default: "hg38", inputBinding: { position: 2, prefix: "--genome_version" } }
  output_basename: { type: string, doc: "String that will be used in the output filenames" }
outputs:
  autopvs1_tsv: { type: File, outputBinding: { glob: '*.autopvs1.tsv'} }