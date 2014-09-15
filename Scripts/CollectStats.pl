#!/usr/bin/perl
use warnings;
use strict;

my $dir = $ARGV[0];

my $samtools = "/Installs/samtools";

opendir(DIR,$dir) || die "Cannot open $dir!\n$!\nexiting...\n";
my @files = grep {/FSR14/} readdir DIR;
closedir DIR;

foreach my $file (@files){
	my $path=$file."/Alignments.sorted.bam";
	my $command = $samtools ." view $path 21119 | wc -l ";
	open(CMD,"-|",$command);
	my @output=<CMD>;
	close CMD;
	chomp $output[0];
	my $hit = $output[0];
	$command = $samtools ." view -F 4 $path | wc -l ";
	open(CMD,"-|",$command);
	@output=<CMD>;
	close CMD;
	chomp $output[0];
	my $total = $output[0];
	my $ratio=$hit/$total;
	print $file.",$total,$hit,$ratio\n";
}
