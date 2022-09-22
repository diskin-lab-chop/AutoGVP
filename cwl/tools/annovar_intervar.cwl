cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-annovar-intervar
label: ANNOVAR for InterVar
doc: |
  "This is a limited application of running annovar with intervar settings"
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
      && perl /annovar/table_annovar.pl

inputs:
  input_vcf: { type: File, secondaryFiles: [.tbi], doc: "VCF file (with associated index) to be annotated", inputBinding: { position: 1} }
  annovar_db: { type: File, doc: "Annovar Database with at minimum required resources to InterVar", inputBinding: { position: 0 }}
  annovar_db_str: { type: string, doc: "Name of dir created when annovar db is un-tarred", inputBinding: { position: 1 }}
  buildver:  { type: ['null', { type: enum, symbols: ["hg38","hg19","hg18"], name: "buildver" } ], doc: "Genome reference build version",
    default: "hg38", inputBinding: { prefix: "--buildver" } }
  remove: { type: 'boolean?', doc: "Remove intermediate files", default: true,
    inputBinding: {position: 1, prefix: "--remove" } }
  output_basename: { type: string, doc: "String that will be used in the output filenames. Be sure to be consistent with this as InterVar will use this too",
    inputBinding: {position: 1, prefix: "--out" } }
  protocol: { type: 'string?', doc: "csv string of databases within `annovar_db` cache to run",
    default: "refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,dbnsfp42a,clinvar_20210501,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene",
    inputBinding: { position: 1, prefix: "--protocol"} }
  operation: { type: 'string?', doc: "csv string of how to treat each listed protocol",
    default: "g,f,f,f,f,f,f,f,r,g,g", inputBinding: { position: 1, prefix: "--operation"} }
  nastring: { type: 'string?', doc: "character used to represent missing values",
    default: '.', inputBinding: { position: 1,  prefix: "--nastring"} }
  otherinfo: { type: 'boolean?', doc: "print out otherinfo (information after fifth column in queryfile)",
    default: true, inputBinding: {position: 1, prefix: "--otherinfo" } }
  threads: { type: 'int?', doc: "Num threads to use to process filter inputs",
    default: 8, inputBinding: { position: 1, prefix: "--threads"} }
outputs:
  annovar_txt: { type: File, outputBinding: { glob: '*._multianno.txt'} }