use warnings;
use strict;

my %info;

open(A, "GENCODE_GRCh37.p12_Ensembl73.txt");
while(<A>){
    chomp;
    my ($ensid,$descr,$gene) = split(/\t/);
    $info{$ensid} = $descr . "\t" . $gene;
}

open(B, "GTEx_Tissue-Wide_local_h2_se.txt");
open(OUT, ">GTEx_Tissue-Wide_local_h2_se_geneinfo.txt");

while(defined(my $line=<B>)){
    my ($e) = split(/\t/,$line);
    my ($ensid) = split(/\./,$e);
    my $geneinfo;
    if(defined($info{$ensid})){
	$geneinfo = $info{$ensid};
    }else{
	$geneinfo = "NA\tNA";
    }
    print OUT "$geneinfo\t$line";
}
