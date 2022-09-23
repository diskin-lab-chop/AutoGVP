cwlVersion: v1.2
class: CommandLineTool
id: kfdrc-intervar
label: InterVar Pathogenicity Classifier
doc: |
  "This tool runs InterVar"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 16000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'migbro/intervar:2.2.1'
  - class: InitialWorkDirRequirement
    listing: [$(inputs.input_ann)]

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      tar --use-compress-program="pigz -p 8" -xf
  - position: 1
    shellQuote: false
    valueFrom: >-
      && tar --use-compress-program="pigz -p 8" -xf
  - position: 2
    shellQuote: false
    valueFrom: >-
      && cp /annovar/*.pl ./
      && python3 /Intervar.py
  - position: 99
    shellQuote: false
    valueFrom: |
      1>&2


inputs:
  input_ann: { type: File, secondaryFiles: ['.tbi?'], doc: "VCF file (with associated index) or ANNOVAR-formatted file to be classified", inputBinding: { position: 2, prefix: "-i"} }
  annovar_db: { type: File, doc: "Annovar Database with at minimum required resources to InterVar", inputBinding: { position: 0 }}
  annovar_db_str: { type: string, doc: "Name of dir created when annovar db is un-tarred", inputBinding: { position: 2, prefix: "-d" }}
  intervar_db: { type: File, doc: "InterVar Database from git repo + mim_genes.txt", inputBinding: { position: 1 }}
  buildver: { type: ['null', { type: enum, symbols: ["hg38","hg19","hg18"], name: "buildver" } ], doc: "Genome reference build version",
    default: "hg38", inputBinding: { position: 2, prefix: "-b" } }
  input_type:  { type: ['null', { type: enum, symbols: ["AVinput","VCF","VCF_m"], name: "input_type" } ], doc: "Input variant type to classify",
    default: "AVinput", inputBinding: { position: 2, prefix: "--input_type" } }
  skip_annovar: { type: boolean, doc: "If giving properly foramtted and annotated ANNOVAR input, skip ANNOVAR run", inputBinding: { position: 2, prefix: "--skip_annovar"}}
  intervar_db_str: { type: string, doc: "Name of dir created when intervar db is un-tarred", inputBinding: { position: 2, prefix: "-t"}}
  output_basename: { type: string, doc: "String that will be used in the output filenames. Be sure to be consistent if a previous run of ANNOVAR was used",
    inputBinding: {position: 2, prefix: "-o" } }

outputs:
   intervar_scored: { type: File, outputBinding: {glob: "*_multianno.txt.intervar"}}