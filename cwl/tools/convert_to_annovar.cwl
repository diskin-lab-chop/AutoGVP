cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-annovar-convert
label: Convert to ANNOVAR FMT
doc: |
  "Convert vcf to annovar input format. Useful for repeat runs of annovar on the same file, especially if it is massive."
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'migbro/intervar:2.2.1'

baseCommand: [perl]
arguments:
  - position: 1
    shellQuote: true
    valueFrom: >-
      /annovar/convert2annovar.pl
      --allsample
      --withfreq
      --format vcf4
      --outfile $(inputs.output_basename).avinput

inputs:
  input_vcf: { type: File, secondaryFiles: [.tbi], doc: "VCF file (with associated index) to be annotated", inputBinding: { position: 1} }
  output_basename: { type: string, doc: "String that will be used in the output filenames" }

outputs:
  vcf_to_gz_avinput: { type: File, outputBinding: { glob: '*.avinput'} }
