#!/usr/bin/perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use Configuration;
use Tools;


my $configFile=$ARGV[0];
my $inDir=$ARGV[1];
my $outDir=$ARGV[2];
my $length=$ARGV[3];

die "usage: perl $0 <QC config file> <input directory (must be *.fastq files)> <output directory> <length of 5\' trim> <threads>\n\n" unless $configFile && $inDir && $outDir && $length;

my $q = Thread::Queue->new();
my $config = Configuration->new($configFile);
my $threads = $config->get("OPTIONS","Threads");
my $script=$config->get("PATHS","FivePrimeTrimmer");

my @files=grep {m/fastq$/} @{Tools->LoadDir($inDir)};
foreach my $file (@files){
	warn "processing $file\n";
	my $Ipath=$inDir."/$file";
	my $Opath=$outDir."/$file";
	$Opath=~s/fastq/5pTrim\.fastq/;
	my $command="perl $script $Ipath $Opath $length";
	$q->enqueue($command);
}
for(my $i=0;$i<=$threads;$i++){
	my $thr=threads->create(\&workerThread);
}
while(threads->list()>0){
	my @thr=threads->list();
	$thr[0]->join();
}

sub workerThread {
	while(my $command=$q->dequeue_nb()){
		warn "$command dequeue'd\n";
		`$command`;
	}
	return 1;
}
