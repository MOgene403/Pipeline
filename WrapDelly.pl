#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/Lib";

use Tools;
use Configuration;

use threads;
use threads::shared;
use Thread::Queue;

my $q = Thread::Queue->new();
die "usage: perl $0 <Config file>\n\n" unless $#ARGV==0;
my $Config = Configuration->new($ARGV[0]);

my $nThreads = $Config->get("OPTIONS","Threads");

warn "Recognizing $nThreads as max threading...\n";

my $ref=$Config->get("PATHS","reference_dir");
warn "Finding Vectors...\n";
my $vecDir = $Config->get("PATHS","vector_dir");
my @LineNo = $Config->getAll("VECTORS");

foreach my $i (@LineNo){
      $q->enqueue($i);
}
for(my$i=0;$i<1;$i++){
      my $thr=threads->create(\&worker);
}
while(threads->list()>0){
      my @thr=threads->list();
      $thr[0]->join();
}


sub worker {
	my $TID=threads->tid() -1 ;
	while(my$j=$q->dequeue_nb()){
		my ($R1,$R2)=split(/\,/,$Config->get("DATA",$j));
		my $base = $Config->get("CELL_LINE",$j);
		my $P1=$Config->get("DIRECTORIES","filtered_dir")."/$base.R1.fastq";
		my $P2=$Config->get("DIRECTORIES","filtered_dir")."/$base.R2.fastq";
		my $inputDir = $Config->get("DIRECTORIES","output_dir")."/".$base;
		my $samtools = $Config->get("PATHS","samtools");
		my $bwaRef=$Config->get("DIRECTORIES","output_dir")."/".$Config->get("VECTORS",$j).".ref.fasta";

		my $bwaRoot=$inputDir."/$base.Alignments";
		my $depths =$inputDir."/$base.ContigDepths.txt";

		my $bwaAln=$bwaRoot.".sorted.bam";

		my $RunDelly=$Config->get("PATHS","RunDelly");
		my $delly = $Config->get("PATHS","delly");

		my  $outputDir = $Config->get("DIRECTORIES","output_dir")."/VCFs";
		mkdir $outputDir unless -e $outputDir;
		warn "making $outputDir\n";
		my $rootLoc = $outputDir."/".$j;
		my $cmd="perl $RunDelly $bwaRef $depths 50 $j ".$bwaAln." ".$ARGV[0];
		warn $cmd."\n";
		`$cmd`;
		#print $cmd."\n";
#usage: perl RunDelly.bestN.pl <reference fasta file> <file of contig IDs to search> <N - number of contigs to search (i.e., '50' for top 50) <id of insertional chromosome> <PE bam><
	}
}

# /home/ec2-user/Store1/bin/delly  -t TRA -o TRA.vcf -q 20 -g TwoChrom.fasta pGC1_Raw.sorted.bam

exit(0);


