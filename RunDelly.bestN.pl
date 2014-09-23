#!/usr/bin/perl
use warnings;
use strict;
use lib '/home/ec2-user/Store1/Tools/Lib';
use FindBin;
use lib "$FindBin::Bin/Lib";
use hdpTools;
use threads;
use threads::shared;
use Thread::Queue;
my $q = Thread::Queue->new();
die "usage: perl $0 <reference fasta file> <file of contig IDs to search> <N - number of contigs to search (i.e., '50' for top 50) <id of insertional chromosome> <PE bam>\n\n" unless $#ARGV==4;

my $file=$ARGV[0];
my @Best=@{hdpTools->LoadFile($ARGV[1])};
my $N  	=$ARGV[2];
my $id	=$ARGV[3];
my $PE	=$ARGV[4];

my @IDs;
foreach my $line (@Best){
	my ($contig,$depth)=split(/\t/,$line);
	push @IDs, $contig;
}

my $temp="./temp";
mkdir $temp unless -e $temp;
my $bkf=$temp."/badkeys.txt";
my $delly="/home/ec2-user/Store1/SVpipe/bin/delly ";


my %Fasta=%{hdpTools->LoadFasta($file)};

my @keys = grep {!/$id/} keys %Fasta;

warn "Avoiding $id!\n";

for(my$i=0;$i<=$N;$i++){
      $q->enqueue($i);
}
for(my$i=0;$i<=7;$i++){
      my $thr=threads->create(\&worker);
}
while(threads->list()>0){
      my @thr=threads->list();
      $thr[0]->join();
}


sub worker {
	my $TID=threads->tid() -1 ;
	while(my$j=$q->dequeue_nb()){
		my $key=$IDs[$j];
		my @badkeys;
		for(my$i=0;$i<=$#keys;$i++){
			push @badkeys, $keys[$i] unless $keys[$i] eq $key;
		}
		my $bk=$bkf.".$j";
		hdpTools->printToFile($bk,\@badkeys);
		my $tfile="Seq".$key;
		$tfile=~s/\|/-/g;
		$tfile.=".vcf";
	#	my $cmd=$delly." -t TRA -o $tfile -q 30 -g $file -x $bk $PE";
		my $cmd=$delly." -t TRA -o $tfile -q 30 -x $bk $PE";
		warn $cmd."\n";
		`$cmd`;
	}
}

# /home/ec2-user/Store1/bin/delly  -t TRA -o TRA.vcf -q 20 -g TwoChrom.fasta pGC1_Raw.sorted.bam

exit(0);


