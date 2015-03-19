#!/usr/bin/perl
use warnings;
use strict;

###pull Framingham eQTL SNP lists at various thresholds: p<0.0001 (entire file), FDR<0.05
###For now, only include hapmap2 SNPs because that's what we selected from the DGN imputation


my $framing = '/group/im-lab/nas40t2/haky/Data/Transcriptome/eQTL-results/Framingham/eqtl-gene-1000g-peer-validate_c1-23_c28.csv.gz';

my $hapmap = '/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/hapmapSnpsCEU.list'; ##from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz

open(OUTALL, ">Framingham_trans-eqtl-gene_p0.0001_hapmapSnpsCEU.csv");
open(OUTFDR, ">Framingham_trans-eqtl-gene_fdr0.05_hapmapSnpsCEU.csv");
open(ALLLIST, ">Framingham_trans-eqtl-gene_p0.0001_hapmapSnpsCEU.list");
open(FDRLIST, ">Framingham_trans-eqtl-gene_fdr0.05_hapmapSnpsCEU.list");


open(HAP, "$hapmap");
my %haplist;
while(<HAP>){
    chomp;
    my ($snp) = split(/\n/);
    $haplist{$snp} = 1;
}

system("zcat $framing > tmp.framing");

open(A, "tmp.framing");

while(defined(my $line = <A>)){
    my @array = split(/,/,$line);
    my $snp = $array[3];
    my $genesymbol = $array[15];
    my $gene;
    if($genesymbol =~ m/\|/){
	($gene) = split(/|/,$genesymbol);
    }else{
	$gene = $genesymbol;
    }
    my $iscis = $array[17];
    if($iscis eq 'Is_Cis'){
        print OUTALL $line;
	print OUTFDR $line;
    }
    elsif(defined($haplist{$snp}) && $gene =~ m/\w+/){
	if($iscis == 0){
	    print OUTALL $line;
	    print ALLLIST "$snp\t$gene\n";
	    my $logfdr = $array[22];
	    my $fdr = 10**$logfdr;
	    if($fdr < 0.05){
		print OUTFDR $line;
		print FDRLIST "$snp\t$gene\n";
	    }
	}
    }   
}

system("rm tmp.framing");
system("gzip Framingham_trans-eqtl-gene*csv");
