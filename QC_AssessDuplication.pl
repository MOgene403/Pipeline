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
#M00645:116:000000000-ABTF3:1:1107:25684:9168    147     9       11501749        0       236M    =       11501606        -379    TATCACAGAACTAAGTGCTTGTGCGGAATCAGTACTGGCTTTTGTTTGGTGGAGGATCAATACTTGCTTTTGTTTGGGGGTGGCAACTGTTTTGCTATAAGATTCCATGTGTTCCTGTTGAGATGAATCATATATAGTATAGCTGCATACTACAAATCTGTTTTTCAAATTTAGGTTGCTTTGGCATGATCTATTTTTTTGTCAGACAGACTTTCTAAGTGGTAGCTCTTGATTTC     FFFFFGEGFFGFGGFGGFBGGGFHGFFHHBFGCFHHHFHFHHHHHHHHHGHHHGGHHHHHG=HGHHHHCHHHHGGGGGHHGHFFHFHHHHHHHFGHHHHHHHHHHHHHHHHHHHHHHHHGHGHHHHHHHHFFFGGHHGHHHHHHHHHHFFHFHHHHHHHGHHHGHHHHHHHHHHHGHHGHHHHHHHHHHHHHHGGHHHHGHHHHHHGCHHHFHHHHHHHHHHHHHHHGHGGGGGGG      NM:i:0  MD:Z:236        AS:i:236        XS:i:236
#M00645:116:000000000-ABTF3:1:1107:24994:9242    147     9       11501749        0       236M    =       11501606        -379    TATCACAGAACTAAGTGCTTGTGCGGAATCAGTACTGGCTTTTGTTTGGTGGAGGATCAATACTTGCTTTTGTTTGGGGGTGGCAACTGTTTTGCTATAAGATTCCATGTGTTCCTGTTGAGATGAATCATATATAGTATAGCTGCATACTACAAATCTGTTTTTCAAATTTAGGTTGCTTTGGCATGATCTATTTTTTTGTCAGACAGACTTTCTAAGTGGTAGCTCTTGATTTC     FFFFBGGFFBFFGGGGGGGGFCDFHEGGHHFHHCHHHHG.HHHHHAHGHHFFGEGBGHHHFFHHFHHGCGHHHGCGGCHHFHFFGBHHGHHHHHHHHHHHGHHHHHHHGHHFHGGHHFHHGHFHHHHHHHGHGHHHHHHGHHHHGGHGHHHHFHHHBHGGHGHHHHHHGHHHHHFHHFFFHFHHFHHHHFHFHGGHFHHHHHHGHHGHCHHFHHHHHHHHHGHHHGFHHGFEFCGF      NM:i:0  MD:Z:236        AS:i:236        XS:i:236
#M00645:116:000000000-ABTF3:1:1107:21080:11822   147     9       11501749        0       236M    =       11501606        -379    TATCACAGAACTAAGTGCTTGTGCGGAATCAGTACTGGCTTTTGTTTGGTGGAGGATCAATACTTGCTTTTGTTTGGGGGTGGCAACTGTTTTGCTATAAGATTCCATGTGTTCCTGTTGAGATGAATCATATATAGTATAGCTGCATACTACAAATCTGTTTTTCAAATTTAGGTTGCTTTGGCATGATCTATTTTTTTGTCAGACAGACTTTCTAAGTGGTAGCTCTTGATTTC^C
#[ec2-user@ip-172-31-39-57 MZSM141720A037A_S20]$ /Installs/samtools view MZSM141720A037A_S20.Alignments.sorted.bam ^C
#
		my $inputDir = $Config->get("DIRECTORIES","output_dir")."/".$base;
		my $samtools = $Config->get("PATHS","samtools");
		my $bwaRef=$Config->get("DIRECTORIES","output_dir")."/".$Config->get("VECTORS",$j).".ref.fasta";

		my $bwaRoot=$inputDir."/$base.Alignments";
		my $depths =$inputDir."/$base.ContigDepths.txt";

		my $bwaAln=$bwaRoot.".sorted.bam";
		my $command = $samtools." view $bwaAln";
		open(CMD,"-|",$command) || die "cannot run $command\n";
		my @output = <CMD>;
		close CMD;
		warn $command."\n";
		my %A;
		foreach my $line (@output){
			my @line = split(/\t/,$line);
			my ($read,$flag,$chr,$pos,@other)=@line;
			if(defined($A{$chr})){
				$A{$chr}{$pos}+=1;
			}else{
				$A{$chr}={};
				$A{$chr}{$pos}+=1;
			}
		}
		my $rep=0;
		my $sin=0;
		foreach my $chr (keys %A){
			foreach my $pos (keys %{$A{$chr}}){
				if($A{$chr}{$pos}>1){
					$rep++;
				}else{
					$sin++;
				}
			}
		}
		warn "for library $base, $rep repeated hits, $sin single hits\n";
	}
}


exit(0);


