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
$nThreads = 3 if $nThreads > 3;

warn "Recognizing $nThreads as max threading...\n";

my $ref=$Config->get("PATHS","reference_dir");
warn "Finding Vectors...\n";
my $vecDir = $Config->get("PATHS","vector_dir");
my @LineNo = $Config->getAll("VECTORS");

foreach my $i (@LineNo){
      $q->enqueue($i);
}
for(my$i=0;$i<$nThreads;$i++){
      my $thr=threads->create(\&worker);
}
while(threads->list()>0){
      my @thr=threads->list();
      $thr[0]->join();
}


sub worker {
	my $TID=threads->tid() -1 ;
	while(my$j=$q->dequeue_nb()){
		my $path=$Config->get("PATHS","vector_dir")."/".$Config->get("VECTORS",$j);
		my $reference = $Config->get("PATHS","reference_dir")."/".$Config->get("BACKGROUNDS",$j);
		warn "Loading $path ...\n";
		my $outputDir = $Config->get("PATHS","temp_dir");
		mkdir $outputDir unless -e $outputDir;
		my $file=$outputDir."/".$j.".ref.fasta";
		my $cmd="cat $reference $path > $file";
		warn $cmd."\n";	
		`$cmd`;
		$cmd=$Config->get("PATHS","bwa")." index $file";
		warn $cmd."\n";
		`$cmd`;
	}
}

# /home/ec2-user/Store1/bin/delly  -t TRA -o TRA.vcf -q 20 -g TwoChrom.fasta pGC1_Raw.sorted.bam

exit(0);


