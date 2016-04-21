#!/usr/bin/perl
use strict;
use warnings;
use Compress::Zlib;

my $geudir = "/group/im-lab/nas40t2/Data/Transcriptome/GEUVADIS/dosagefiles/";


for(my $i = 1; $i <= 22; $i++){
    my $dosefile = ${geudir} . "chr" . $i . ".dosage.gz";
    my $gzdose = gzopen($dosefile, "rb") or die "Cannot open $dosefile: $gzerrno\n";

    open(OUT, ">geu_dosage_1KGP_MAF0.01/chr" . $i . ".dosage");

    my $line;

    while($gzdose->gzreadline($line)){
	my ($chr, $snp, $bp, $af2, $a1, $a2, @rest) = split(/\t/, $line);
	if($af2 > 0.01 && $af2 < 0.99){
	    print OUT $line;
	}
    }	

    $gzdose->gzclose();
    system("gzip geu_dosage_1KGP_MAF0.01/chr" . $i . ".dosage");
}
