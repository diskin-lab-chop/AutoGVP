#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27 19:13
# Modified: 2022/10/05
# Mod_auther: Miguel Brown

import os
import re
import sys
import random
import string
import argparse
from collections import namedtuple
import pysam
from read_data_mod import trans_gene, gene_trans, gene_alias
from read_data_mod import transcripts_hg38
from pvs1 import PVS1
from utils import get_transcript, vep_consequence_trans, VCFRecord
__version__ = 'v0.2.0'


lof_type = ['frameshift', 'nonsense', 'splice-5', 'splice-3', 'init-loss']
vep_lof_list = ['frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 'start_lost']
VAR = namedtuple('VAR', ('varid', 'gene', 'trans', 'canonical', 'pick', 'record'))


class AutoPVS1:
    """Run AutoPVS1"""
    @staticmethod
    def get_vcfrecord(vcf):
        if all(hasattr(vcf, attr) for attr in ['chrom', 'pos', 'ref', 'alt']):
            return vcf
        elif re.match(r'^((?:chr)?[\dXYMT]{1,2})-(\d+)-([ATCG]+)-([ATCG]+)$', vcf, re.I):
            v = vcf.split("-")
            v[0] = re.sub('chr', '', v[0], flags=re.I)
            return VCFRecord(v[0], int(v[1]), v[2], v[3])
        else:
            raise TypeError("Wrong VCF Record")


    def run_pvs1(self):
        pvs1 = PVS1(self.vcfrecord, self.consequence, self.hgvs_c, 
                    self.hgvs_p, self.transcript, self.genome_version)
        return pvs1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vep_vcf', help='VCF classify')
    parser.add_argument('--genome_version', help='Currently only hg38 or GRCh38 supported')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()

    genome_version = args.genome_version
    in_vcf = pysam.VariantFile(args.vep_vcf, threads=8)

    # Use VEP PICK field to choose a representative transcript
    csq_fields = in_vcf.header.info['CSQ'].description.replace("Consequence annotations from Ensembl VEP. Format: ", "").split("|")

    # get index of PICK field
    p_idx = csq_fields.index("PICK")

    print ("vcf_id",'SYMBOL','Feature','trans_name','consequence', 'strength_raw', 'strength','criterion', sep="\t")

    for record in in_vcf.fetch():
        # Parse CSQ for PICKed transcript
        picked_csqs = []
        try:
            for csq in record.info['CSQ']:
                info_list = csq.split("|")
                if info_list[p_idx] == '1':
                    picked_csqs.append(csq)
        except Exception as e:
            sys.stderr.write(str(e) + '\n')
            sys.stderr.write(record.contig + "\t" + str(record.pos) + '\n')
            continue
        if len(picked_csqs) > 1:
            sys.stderr.write("WARN: More than 1 PICK detected\n" + "\n".join(picked_csqs) + "\n")
        # Populate dict with CSQ key-value pairs
        info = {}
        first_picked = picked_csqs[0].split('|')
        for i in range(len(first_picked)):
            info[csq_fields[i]] = first_picked[i]
            # adjust this filed to match behavior of original auto pvs1 script
            if csq_fields[i] == 'HGVSp':
                info[csq_fields[i]] = first_picked[i].replace('%3D', '=')

        vcfrecord = VCFRecord(record.contig.replace('chr', ''), str(record.pos), record.ref, record.alts[0])
        transcript = get_transcript(info['Feature'], transcripts_hg38)

        consequence = vep_consequence_trans(info['Consequence'])
        vcf_id = "-".join([vcfrecord.chrom, str(vcfrecord.pos), vcfrecord.ref, vcfrecord.alt])
        if consequence in lof_type and transcript:
            lof_pvs1 = PVS1(vcfrecord, consequence, info['HGVSc'], info['HGVSp'], transcript, genome_version)
            trans_name = lof_pvs1.transcript.full_name
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  trans_name,
                  lof_pvs1.consequence,
                  lof_pvs1.strength_raw.name,
                  lof_pvs1.strength.name,
                  lof_pvs1.criterion,
                  sep="\t")
        elif consequence in lof_type:
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'not_canonical',
                  consequence,
                  'Unmet',
                  'Unmet',
                  'na',
                  sep="\t")
        else:
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'not_lof',
                  consequence,
                  'Unmet',
                  'Unmet',
                  'na',
                  sep="\t")

        '''
        if info['Function'] in ['splice-5', 'splice-3']:
            pvs1 = PVS1(vcfrecord, info['Function'], info['pHGVS1'], info['MIMInheritance'], transcript)
            splice = Splicing(vcfrecord, transcript)
            print(vcf_id,
                  pvs1.function,
                  pvs1.strength_raw.name,
                  pvs1.strength.name,
                  splice.is_undergo_NMD,
                  splice.has_cryptic_splice_site,
                  splice.is_exon_skipping,
                  splice.preserves_reading_frame,
                  splice.is_critical_to_protein_func,
                  splice.variant_removes_10_percent_of_protein,
                  splice.is_critical_to_protein_func_detail,
                  sep="\t")
        '''


if __name__ == '__main__':
    main()
