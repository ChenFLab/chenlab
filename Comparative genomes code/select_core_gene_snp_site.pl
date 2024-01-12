#!/usr/bin/perl
use strict;
use warnings;

open IN,"<$ARGV[0]" or die $!; #输入文件step_4_core_snp/../roary_core_gene_site.txt
open IN1,"<$ARGV[1]" or die $!; #输入文件2，step3/4_all/snp.site.uniq.indel.filter.txt
open OUT,">$ARGV[2]" or die $!; #./rm_no_core_snp.txt
open OUT1,">$ARGV[3]" or die $!; #total_core_snp.txt
my %core;
while (<IN>) {
	chomp $_;
	my @lines=split /\t/,$_;
	for (my $i=$lines[0];$i<=$lines[1];$i++){
		$core{$i}=1;
	}
}

my %hash;
my @array;

while (<IN1>) {
	chomp $_;
	
	@array=split /\t/,$_;
	if (exists $core{$array[0]}) {
		print OUT1 "$_\n";
	}
	else {
		print OUT "$_\n";
	}	
}

close IN;
close IN1;
close OUT;
close OUT1;