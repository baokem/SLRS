#!/usr/bin/perl
use strict;
unless (@ARGV) {die "$0 <fastq1> <fastq2>\n";}
my $fq1=shift;
my $fq2=shift;


open (FQ1,"$fq1") or die ("can not open $fq1");
open (FQ2,"$fq2") or die ("can not open $fq2");

while (my $line1_1=<FQ1>){
	my $line1_2=<FQ1>;
	my $line1_3=<FQ1>;
	my $line1_4=<FQ1>;
	
	my $line2_1=<FQ2>;
	my $line2_2=<FQ2>;
	my $line2_3=<FQ2>;
	my $line2_4=<FQ2>;
	print $line1_1,$line1_2,$line1_3,$line1_4,$line2_1,$line2_2,$line2_3,$line2_4;
}
close FQ1;close FQ2;
