cwlVersion: v1.2
class: Workflow
id: intervar-classification-wf
label: InterVar Classification Workflow
doc: |
  Takes in a VCF file, and converts to ANNOVAR table format, annotates with supporting data needed for classification using ANNOVAR, then classifies variants

inputs:
  input_vcf: { type: File, secondaryFiles: [.tbi], doc: "VCF file (with associated index) to be annotated" }
  annovar_db: { type: File, doc: "Annovar Database with at minimum required resources to InterVar" }
  annovar_db_str: { type: string, doc: "Name of dir created when annovar db is un-tarred" }
  buildver:  { type: ['null', { type: enum, symbols: ["hg38","hg19","hg18"], name: "buildver" } ], doc: "Genome reference build version",
    default: "hg38" }
  output_basename: { type: string, doc: "String that will be used in the output filenames. Be sure to be consistent with this as InterVar will use this too" }
  annovar_protocol: { type: 'string?', doc: "csv string of databases within `annovar_db` cache to run",
    default: "refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,dbnsfp42a,clinvar_20210501,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene" }
  annovar_operation: { type: 'string?', doc: "csv string of how to treat each listed protocol",
    default: "g,f,f,f,f,f,f,f,r,g,g" }
  annovar_nastring: { type: 'string?', doc: "character used to represent missing values", default: '.' }
  annovar_otherinfo: { type: 'boolean?', doc: "print out otherinfo (information after fifth column in queryfile)", default: true }
  annovar_threads: { type: 'int?', doc: "Num threads to use to process filter inputs", default: 8 }
  annovar_vcfinput: { type: 'boolean?', doc: "Annotate vcf and generate output file as vcf", default: false }
  intervar_db: { type: File, doc: "InterVar Database from git repo + mim_genes.txt" }
  intervar_db_str: { type: string, doc: "Name of dir created when intervar db is un-tarred" }
outputs:
  intervar_classification: { type: File, outputSource: intervar_classify/intervar_scored}
  annovar_vcfoutput: { type: 'File?', outputSource: run_annovar/vcf_output}

steps:
  convert_to_annovar:
    run: ../tools/convert_to_annovar.cwl
    in:
      input_vcf: input_vcf
      output_basename: output_basename
    out: [vcf_to_avinput]
  run_annovar:
    run: ../tools/annovar_intervar.cwl
    in:
      av_input: convert_to_annovar/vcf_to_avinput
      annovar_db: annovar_db
      annovar_db_str: annovar_db_str
      buildver: buildver
      output_basename: output_basename
      protocol: annovar_protocol
      operation: annovar_operation
      nastring: annovar_nastring
      otherinfo: annovar_otherinfo
      threads: annovar_threads
      vcfinput: annovar_vcfinput
    out: [annovar_txt, vcf_output]
  intervar_classify:
    run: ../tools/intervar.cwl
    in:
      input_ann: run_annovar/annovar_txt
      annovar_db: annovar_db
      annovar_db_str: annovar_db_str
      intervar_db: intervar_db
      buildver: buildver
      output_basename: output_basename
      intervar_db_str: intervar_db_str
      skip_annovar:
        valueFrom: "${return true;}"
    out: [intervar_scored]



