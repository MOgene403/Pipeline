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

die "usage: perl $0 <QC config file> <input directory (must be *.fastq files)> <output directory> <threads>\n\n" unless $configFile && $inDir && $outDir;

my $config=Configuration->new($configFile);

my $threads = $config->get("OPTIONS","managerThreads") * $config->get("OPTIONS","Threads");

my $q = Thread::Queue->new();

my $script=$config->get("PATHS","fastq_quality_trimmer");
my $opts="-t 30 -Q 33";


my @files=grep {m/fastq$/} @{Tools->LoadDir($inDir)};
unless(defined($files[0])){
	@files=grep {m/fq$/} @{Tools->LoadDir($inDir)};
}

foreach my $file (@files){
	my $Ipath=$inDir."/$file";
	my $Opath=$outDir."/$file";
	$Opath=~s/fastq/3pQtrim\.fastq/;
	$Opath=~s/fq/3pQtrim\.fastq/;
	my $command="$script $opts -i $Ipath -o $Opath";
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
