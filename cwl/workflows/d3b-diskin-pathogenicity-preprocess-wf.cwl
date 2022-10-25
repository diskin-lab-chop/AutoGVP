cwlVersion: v1.2
class: Workflow
id: d3b-diskin-pathogenicity-preprocess-wf
label: Pathogenicity Preprocessing Workflow
doc: |
  Runs InterVar and autoPVS1 on a VEP-annotated VCF
requirements:
- class: SubworkflowFeatureRequirement

inputs:
  # Common
  vep_vcf: { type: File, secondaryFiles: [.tbi], doc: "VCF file (with associated index) to be annotated" }
  buildver:  { type: ['null', { type: enum, symbols: ["hg38","hg19","hg18"], name: "buildver" } ], doc: "Genome reference build version",
    default: "hg38" }
  output_basename: { type: string, doc: "String that will be used in the output filenames. Be sure to be consistent with this as InterVar will use this too" }
  # InterVar
  annovar_db: { type: File, doc: "Annovar Database with at minimum required resources to InterVar" }
  annovar_db_str: { type: 'string?', doc: "Name of dir created when annovar db is un-tarred", default: "annovar_humandb_hg38_intervar" }
  annovar_protocol: { type: 'string?', doc: "csv string of databases within `annovar_db` cache to run",
    default: "refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,dbnsfp42a,clinvar_20210501,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene" }
  annovar_operation: { type: 'string?', doc: "csv string of how to treat each listed protocol",
    default: "g,f,f,f,f,f,f,f,r,g,g" }
  annovar_nastring: { type: 'string?', doc: "character used to represent missing values", default: '.' }
  annovar_otherinfo: { type: 'boolean?', doc: "print out otherinfo (information after fifth column in queryfile)", default: true }
  annovar_threads: { type: 'int?', doc: "Num threads to use to process filter inputs", default: 8 }
  intervar_db: { type: File, doc: "InterVar Database from git repo + mim_genes.txt" }
  intervar_db_str: { type: 'string?', doc: "Name of dir created when intervar db is un-tarred", default: "intervardb" }
  # autoPVS1
  autopvs1_db: { type: File, doc: "git repo files plus a user-provided fasta reference"}
  autopvs1_db_str: { type: 'string?', doc: "Name of dir created when annovar db is un-tarred", default: "data"}

outputs:
  intervar_classification: { type: File, outputSource: run_intervar/intervar_classification}
  autopvs1_tsv: { type: File, outputSource: run_autopvs1/autopvs1_tsv }

steps:
  run_intervar:
    run: intervar_classificatiion_wf.cwl
    in:
        input_vcf: vep_vcf
        annovar_db: annovar_db
        annovar_db_str: annovar_db_str
        buildver: buildver
        output_basename: output_basename
        annovar_protocol: annovar_protocol
        annovar_operation: annovar_operation
        annovar_nastring: annovar_nastring
        annovar_otherinfo: annovar_otherinfo
        annovar_threads: annovar_threads
        intervar_db: intervar_db
        intervar_db_str: intervar_db_str
    out: [intervar_classification]

  run_autopvs1:
    run: ../tools/autopvs1.cwl
    in:
        autopvs1_db: autopvs1_db
        autopvs1_db_str: autopvs1_db_str
        vep_vcf: vep_vcf
        genome_version:  buildver
        output_basename: output_basename
    out: [autopvs1_tsv]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: sbg:maxNumberOfParallelInstances
  value: 2
