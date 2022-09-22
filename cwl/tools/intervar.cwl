cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-intervar
label: ANNOVAR for InterVar
doc: |
  "This tool runs InterVar"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'migbro/intervar:2.2.1'

baseCommand: []
arguments:
  - position: 0
    shellQuote: true
    valueFrom: >-
      tar --use-compress-program="pigz -p $(inputs.threads)" -xf
  - position: 1
    shellQuote: true
    valueFrom: >-
      && tar --use-compress-program="pigz -p $(inputs.threads)" -xf
  - position: 2
    shellQuote: true
    valueFrom: >-
      && cp /annovar/*.pl ./
      && python3 /Intervar.py -b hg38 -i intervar_Test2.hg38_multianno.txt --input_type=AVinput --database_intervar /WORK/tools/InterVar/intervardb --database_locat intervar_humandb_hg38 --skip_annovar -o intervar_Test2

inputs:
  input_vcf: { type: File, secondaryFiles: [.tbi], doc: "VCF file (with associated index) to be annotated", inputBinding: { position: 1, prefix: "-i"} }
  annovar_db: { type: File, doc: "Annovar Database with at minimum required resources to InterVar", inputBinding: { position: 0 }}
  annovar_db_str: { type: string, doc: "Name of dir created when annovar db is un-tarred", inputBinding: { position: 2, prefix: "-d" }}
  intervar_db: { type: File, doc: "InterVar Database from git repo + mim_genes.txt", inputBinding: { position: 1 }}
  buildver: { type: ['null', { type: enum, symbols: ["hg38","hg19","hg18"], name: "buildver" } ], doc: "Genome reference build version",
    default: "hg38", inputBinding: { position: 2, prefix: "-b" } }
  input_type:  { type: ['null', { type: enum, symbols: ["AVinput","VCF","VCF_m"], name: "input_type" } ], doc: "Input variant type to classify",
    default: "AVinput", inputBinding: { position: 2, prefix: "--input_type" } }
  skip_annovar: { type: boolean, doc: "If giving properly foramtted and annotated ANNOVAR input, skip ANNOVAR run", inputBinding: { position: 2, prefix: "--skip_annovar"}}
  intervar_db_str: { type: string, doc: "Name of dir created when intervar db is un-tarred", inputBinding: { position: 2, prefix: "--database_intervar"}}
  output_basename: { type: string, doc: "String that will be used in the output filenames. Be sure to be consistent if a previous run of ANNOVAR was used",
    inputBinding: {position: 2, prefix: "-o" } }

outputs:
   intervar_scored: { type: File, outputBinding: {glob: "*_multianno.txt.intervar"}}