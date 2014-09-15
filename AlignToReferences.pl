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
		my ($R1,$R2)=split(/\,/,$Config->get("GROUPS",$j));
		my $P1=$Config->get("PATHS","output_dir")."/".$j.".R1.fastq";
		my $P2=$Config->get("PATHS","output_dir")."/".$j.".R2.fastq";
		my $outputDir = $Config->get("PATHS","output_dir");
		my $samtools = $Config->get("PATHS","samtools");
		my $bwaRef= $Config->get("PATHS","temp_dir")."/$j.ref.fasta";
		my $bwaRoot=$outputDir."/$j.Alignments";
		my $bwaAln=$bwaRoot.".bam";
		my $cmd=$Config->get("PATHS","bwa")." mem -t $nThreads $bwaRef $P1 $P2 | $samtools view -bS - > $bwaAln";
		warn $cmd."\n";
		`$cmd`;
		my $sorted=$bwaRoot.".sorted";
		$cmd = $samtools." sort $bwaAln $sorted";
		`$cmd`;
		$cmd = $samtools." index ".$sorted.".bam";
		`$cmd`;
		my $depthscript = $Config->get("PATHS","depthScript");
		my $depthout	= $outputDir."/ContigDepths.txt";
		$cmd = $samtools." depth ".$sorted.".bam | perl $depthscript > $outputDir/$j.ContigDepths.txt";
		`$cmd`;
	}
}

# /home/ec2-user/Store1/bin/delly  -t TRA -o TRA.vcf -q 20 -g TwoChrom.fasta pGC1_Raw.sorted.bam

exit(0);


