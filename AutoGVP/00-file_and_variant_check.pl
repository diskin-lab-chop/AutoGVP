#!/usr/bin/perl

use strict;

################################################################################
# 00-file_and_variant_check.pl
# written by Ammar Naqvi
#
# This script checks if required files exist and if the variants match up
#
# usage: perl 00-file_and_variant_check.pl <optional input dir>
################################################################################

my $dir_to_check = "";
my $vcf_file_check = "";
my $intervar_file_check = "";
my $multianno_file_check = "";
my $autopvs1_file_check = "";
my $dir_to_test = $ARGV[0];

if($dir_to_test=~/\w+/){
  $dir_to_check = $dir_to_test;
}
else{
  $dir_to_check = "input";
}

## check to see if input files all exist
if (glob($dir_to_check."/*.vcf*")) {
  print "vcf file exists...\n";
  $vcf_file_check = $dir_to_check."/*.vcf*";

}
else {
  die("FAIL: vcf file does not exist\n");
}

if (glob($dir_to_check."/*.intervar*")) {
  print "intervar file exists...\n";
  $intervar_file_check = $dir_to_check."/*.intervar*";
}
else {
  die("FAIL: intervar file does not exist\n");
}

if (glob($dir_to_check."/*hg38_multianno*")) {
  print "multianno file exists...\n";
  $multianno_file_check = $dir_to_check."/*hg38_multianno*";
}
else {
  die("FAIL: multianno file does not exist\n");
}

if (glob($dir_to_check."/*autopvs1.txt")) {
  print "autopvs1 file exists...\n";
  $autopvs1_file_check = $dir_to_check."/*autopvs1*";

}
else {
  die("FAIL: autopvs1 file does not exist\n");
}


if (glob($dir_to_check."/*variant_summary*")) {
  print "variant_summary file exists...\n";
}
else {
  die("FAIL: variant_summary file does not exist\n");
}

if (glob($dir_to_check."/*submission_summary*")) {
  print "submission_summary file exists...\n";
}
else {
  die("FAIL: submission summary file does not exist\n");
}

## check to see variants match up with intervar, autopvs1 and vcf files
my $vcf_file = <\\$vcf_file_check>;
my $intervar_file = <\\$vcf_file_check>;
my $multianno_file = <\\$multianno_file_check>;
my $autopvs1_file = <\\$autopvs1_file_check>;

my %vcf_variants;
my %multianno_variants;
my %autopvs1_variants;

if ($vcf_file=~/gz/){
  open(FIL, "gunzip -c  $vcf_file |") || die("Cannot Open File $vcf_file")
}
else{
  open(FIL,$vcf_file) || die("Cannot Open File $vcf_file");
}
while(<FIL>)
{
  chomp;
  next if $_=~/^#/;
  my @cols = split;
  my $chr = $cols[0];
  my $start = $cols[1];
  my $loc = $chr."-".$start;
  $vcf_variants{$loc} = $loc;

}
close(FIL);

if ($multianno_file=~/gz/){
  open(FIL, "gunzip -c  $multianno_file |") || die("Cannot Open File $multianno_file")
}
else{
  open(FIL,$multianno_file) || die("Cannot Open File $multianno_file");
}


while(<FIL>)
{
  chomp;
  next if $_=~/^#/;
  next if $_=~/Start/;

  if($_=~/(chr(\d+)|[XY])\t(\d+)\t(rs\d+|\.)\t/)
  {
    my @cols = split;
    my $chr = $1;
    my $start = $2;
    my $loc = $chr."-".$start;
    $multianno_variants{$loc} = $loc;
  }
}
close(FIL);

if ($autopvs1_file=~/gz/){
  open(FIL, "gunzip -c  $autopvs1_file |") || die("Cannot Open File $autopvs1_file")
}
else{
  open(FIL,$autopvs1_file) || die("Cannot Open File $autopvs1_file");
}

while(<FIL>)
{
  chomp;
  next if $_=~/^#/;
  next if $_=~/strength_raw/;
  my @cols = split "-";
  my $chr = "chr".$cols[0];
  my $start = $cols[1];
  my $loc = $chr."-".$start;

  $autopvs1_variants{$loc} = $loc;
}
close(FIL);

## match up variants to make sure they are consistent in files
for ( keys %vcf_variants ) {
    unless ( exists $multianno_variants{$_} ) {
        die("$_: not found in multianno");
        next;
    }
    unless ( exists $autopvs1_variants{$_} ) {
        die("$_: not found in autopvs1");
        next;
    }
}
