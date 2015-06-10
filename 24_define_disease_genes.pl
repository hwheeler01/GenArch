use warnings;
use strict;

my $pthresh = $ARGV[1];

open(BD, ">OTD_enrichment/gwas.BD.$pthresh.tsv");
open(CAD, ">OTD_enrichment/gwas.CAD.$pthresh.tsv");
open(HT,">OTD_enrichment/gwas.HT.$pthresh.tsv");
open(T1D,">OTD_enrichment/gwas.T1D.$pthresh.tsv");
open(T2D,">OTD_enrichment/gwas.T2D.$pthresh.tsv");
open(CD,">OTD_enrichment/gwas.CD.$pthresh.tsv");
open(RA,">OTD_enrichment/gwas.RA.$pthresh.tsv");
open(ALL, ">OTD_enrichment/gwas.ALL.$pthresh.tsv"); ##all disease/trait genes in the gwas catalog

my %bd;
my %cad;
my %ht;
my %t1d;
my %t2d;
my %cd;
my %ra;
my %all;

open(FILE, "$ARGV[0]");

while(defined(my $line = <FILE>)){
    my ($a, $b, $c, $d, $e, $f, $g, $trait, $i, $j, $k, $l, $m, $report, $map, $p, $q, $r, $s, $t, $u, $v, $w, $x, $y, $z, $aa, $pval) = split(/\t/,$line);
#    print "$trait\n";
    if($pval =~ m/\d+/ && $pval <= $pthresh){
	$all{$report} = 1;
	if($map =~ m/ - /){
	    my ($g1,$g2) = split(/ - /,$map);
	    $all{$g1} = 1;
	    $all{$g2} = 1;
	}else{
	    $all{$map} = 1;
	}
	if($trait =~ m/Bipolar/){
	    $bd{$report} = 1;
	    if($map =~ m/ - /){
		my ($g1,$g2) = split(/ - /,$map);
		$bd{$g1} = 1;
		$bd{$g2} = 1;
	    }else{
		$bd{$map} = 1;
	    }
	}
	elsif($trait =~ m/Coronary/){
	    $cad{$report} = 1;
	    if($map =~ m/ - /){
		my ($g1,$g2) = split(/ - /,$map);
		$cad{$g1} = 1;
		$cad{$g2} = 1;
	    }else{
		$cad{$map} = 1;
	    }
	}
	elsif($trait =~ m/Hypertension/){
	    $ht{$report} = 1;
	    if($map =~ m/ - /){
		my ($g1,$g2) = split(/ - /,$map);
            $ht{$g1} = 1;
		$ht{$g2} = 1;
	    }else{
		$ht{$map} = 1;
	    }
	}
	elsif($trait eq "Type 1 diabetes"){
	    $t1d{$report} = 1;
	    if($map =~ m/ - /){
		my ($g1,$g2) = split(/ - /,$map);
		$t1d{$g1} = 1;
            $t1d{$g2} = 1;
	    }else{
		$t1d{$map} = 1;
	    }
	}
	elsif($trait =~ m/Type 2 diabetes/ && $trait ne "Type 2 diabetes nephropathy"){
	    $t2d{$report} = 1;
	    if($map =~ m/ - /){
		my ($g1,$g2) = split(/ - /,$map);
		$t2d{$g1} = 1;
		$t2d{$g2} = 1;
	    }else{
            $t2d{$map} = 1;
	    }
	}
	elsif($trait =~ m/Crohn/ || $trait =~ m/Ulcerative/){
	    $cd{$report} = 1;
	    if($map =~ m/ - /){
		my ($g1,$g2) = split(/ - /,$map);
		$cd{$g1} = 1;
		$cd{$g2} = 1;
	    }else{
		$cd{$map} = 1;
	    }
	}
	elsif($trait =~ m/Rheumatoid arthritis/){
	    $ra{$report} = 1;
	    if($map =~ m/ - /){
		my ($g1,$g2) = split(/ - /,$map);
		$ra{$g1} = 1;
		$ra{$g2} = 1;
	    }else{
		$ra{$map} = 1;
	    }
	}
    }
}

foreach my $g (sort(keys %bd)){
    if($g eq ""){
	next;
    }
    print BD "$g\n";
}

foreach my $g (sort(keys %cad)){
    if($g eq ""){
        next;
    }
    print CAD "$g\n";
}

foreach my $g (sort(keys %ht)){
    if($g eq ""){
        next;
    }
    print HT "$g\n";
}

foreach my $g (sort(keys %t1d)){
    if($g eq ""){
        next;
    }
    print T1D "$g\n";
}

foreach my $g (sort(keys %t2d)){
    if($g eq ""){
        next;
    }
    print T2D "$g\n";
}

foreach my $g (sort(keys %cd)){
    if($g eq ""){
        next;
    }
    print CD "$g\n";
}

foreach my $g (sort(keys %ra)){
    if($g eq ""){
        next;
    }
    print RA "$g\n";
}

foreach my $g (sort(keys %all)){
    if($g eq ""){
        next;
    }
    print ALL "$g\n";
}
