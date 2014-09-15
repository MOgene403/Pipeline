#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../Lib";
use Configuration;
use Tools;

my $configFile = $ARGV[0];
my $dir	=$ARGV[1];
my $odir=$ARGV[2];
my $config= Configuration->new($configFile);
die "usage: perl $0 <sequence QC configuration file> <directory containing .fastq files> <output dir>\n\n" unless $#ARGV==2;

opendir(DIR,$dir) || die "cannot open directory $dir!\n$!\nexiting...\n";
my @files = grep {m/fastq$/} readdir(DIR);
unless(defined($files[0])){
	seekdir(DIR,0);
	@files=grep {m/fq$/} readdir (DIR);
}
close DIR;

my $cutAdapt= $config->get('PATHS','CutAdapt');

my $command=$cutAdapt." -e 0.1 -n 5 -O 10 --match-read-wildcards -b  GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -b CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b ATCTCGTATGCCGTCTTCTGCTTG -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA -b CAAGCAGAAGACGGCATACGAGAT -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGC -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCACC";

foreach my $file (@files){
	my $path=$dir."/".$file;
	my $out=$file;
	if($out=~m/fastq/){
		$out=~s/fastq/cutadapt\.fastq/;
	}elsif($out=~m/fq$/){
		$out=~s/fq/cutadapt\.fastq/;
	}
	$out=$odir."/".$out;
	my $cmd=$command." $path -o $out";
	warn $cmd."\n";
	`$cmd`;
}
