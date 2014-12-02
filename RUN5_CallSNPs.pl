#!/usr/bin/perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use FindBin;
use lib "$FindBin::Bin/Lib";
use Configuration;
use Tools;

my $configFile=$ARGV[0];

die "usage : perl $0 <config file governing all alignment>\n\n" unless $#ARGV==0;

my $q = Thread::Queue->new();
my $config = Configuration->new($configFile);
my $threads = $config->get("OPTIONS","Threads");
my @Groups = $config->getAll("DATA");
my $mchScriptDir = $FindBin::Bin."/Lib/SysCall_1_1/";
my $mchScript=$mchScriptDir."/SysCall.pl";
for(my $i=0;$i<=$#Groups;$i++){
	warn "enqueuing $i ($Groups[$i])\n";
	$q->enqueue($Groups[$i]);
}

for(my$i=0;$i<1;$i++){
	my $thr=threads->create(\&workerThread);
}
while(threads->list()>0){
	my @thr=threads->list();
	$thr[0]->join();
}


sub workerThread{
	while(my $work=$q->dequeue_nb()){
		my $grp		= $work;
		my $base 	= $config->get("CELL_LINE",$grp);
		my $p1		= $base.".R1.fastq";
		my $p2		= $base.".R2.fastq";
		my $DataDir 	= $config->get("DIRECTORIES","filtered_dir");
		my $OutDir  	= $config->get("DIRECTORIES","output_dir");
		my $RefDir	= $config->get("DIRECTORIES","reference");
		my $TempDir	= $config->get("DIRECTORIES","temp_dir");
		checkDir($DataDir,$OutDir,$RefDir,$TempDir);
		
		my $outputDir 	= $config->get("DIRECTORIES","output_dir")."/".$config->get("CELL_LINE",$grp);
		my $bwaRoot	= $outputDir."/$base.Alignments";
		my $samtools 	= $config->get("PATHS","samtools");
		my $bwa		= $config->get("PATHS","bwa");
		my $bwaRef	= $config->get("DIRECTORIES","output_dir")."/".$config->get("VECTORS",$grp).".ref.fasta";

		my $workThreads = $config->get("OPTIONS","Threads");
		my $bcftools	= $config->get("PATHS","bcftools");
		my $snpRate	= $config->get("OPTIONS","snpRate");
		my $minCov	= $config->get("OPTIONS","minCov");
		my $ins 	= $config->get("INSERTS",$grp);
		my @CurrentSourcePaths;
		my @GarbageCollector;

		my $file1=$DataDir."/".$p1;
		my $file2=$DataDir."/".$p2;
		die "Cannot find read 1 for group: $grp\nFile missing: $file1\nexiting...\n" unless -e $file1;
		die "Cannot find read 2 for group: $grp\nFile missing: $file2\nexiting...\n" unless -e $file2;

		my $index = $OutDir."/".$config->get("VECTORS",$grp).".ref.fasta";

		my $command = "$samtools mpileup -r $ins -F 0.00001 -g -C50 -d 10000000 -f $index $bwaRoot.sorted.bam | $bcftools view -b -m 0.01 -p .99 - | $bcftools view - > $bwaRoot.raw.vcf";
		warn $command."\n";
		`$command`;
		my %H = %{parseResults("$bwaRoot.raw.vcf",$bwaRoot.".filt.vcf",$snpRate,$minCov)};
		collectTheGarbage(@GarbageCollector);
	}
}

sub parseToFinal {
	my $line=shift;
	my @line=split(/\t/,$line);
	my $chr=$line[0];
	my $pos=$line[1];
	my $refBase=$line[3];
	next if $line[4] eq "X";
        my @posAlt =split(/\,/,$line[4]);
	my %I=%{parseInfo($line[7])};
			
}

sub getSysErrors {
	my $file=shift;
	my @file=@{Tools->LoadFile($file.".sys_errors")};
	my $head=shift @file;
	my %H;
	foreach my $line (@file){
		$line=~s/\s.+//;
		my ($chr,$coord)=split(/\:/,$line);
		my $k="$chr-$coord";
		$H{$k}=1;
	}
	return \%H;
}

sub generateMeacham {
	my $file=shift;
	my $out=shift;
	my %H=%{$_[0]};
	my @output;
	my %Fasta=%{Tools->LoadFasta($file)};
	foreach my $pos (keys %H){
		my $chr=$Fasta{$H{$pos}{"Chr"}};
		next if $pos <3;
		my $rpos=$pos-3;  ### -1 for coordinate correct, -2 for start of meacham.
		my $seq=substr($chr,$rpos,5);
		my @seq=split(//,$seq);
		push @output, $H{$pos}{"Chr"}." ".$pos." ".join(" ",@seq);
	}
	return 0 unless scalar(@output)>0;
	Tools->printToFile($out,\@output);
	return 1;
}

sub checkDir {
	for(my$i=0;$i<=$#_;$i++){
		my $dir=$_[$i];
		unless(-e $dir){
			warn "Problem with configuration ($i).\n";
			die "Cannot find directory: $dir\n";
		}
	}
	return 1;
}

sub parseResults {
	my $file=shift;
	my $output=shift;
	my $rate=shift;
	my $minCov=shift;
	my @file=@{Tools->LoadFile($file)};
	my @out;
	my %R;
	foreach my $line (@file){
		my @line=split(/\t/,$line);
		my $chr=$line[0];
		my $pos=$line[1];
		my $refBase=$line[3];
		next if $line=~m/^\#/;
		next if $line[4] eq "X";
		next if $line =~m/INDEL/;
		my @posAlt =split(/\,/,$line[4]);
		my %I=%{parseInfo($line[7])};
		unless(defined $I{"DP"}){
			warn $line."\n";
			die "malformed output from SNP calling!\n";
		}
		unless(defined $I{"I16"}){
			warn $line."\n";
			die "malformed output from SNP calling!\n"
		}
		unless( defined $I{"QS"}){
			warn $line."\n";
			die "malformed output from SNP calling!\n";
		}
		my %r;
		$r{"Chr"}=$chr;
		$r{"Pos"}=$pos;
		$r{"DP"}=$I{"DP"};
		$r{"I16"}=$I{"I16"};
		$r{"QS"}=$I{"QS"};
		$r{"L"}=$line;
		if(checkInfo(\%r,$rate,$minCov)){
			$R{$pos}=\%r;
			push @out, $line;
		}
	}
#	Tools->printToFile($output,\@out);
	return \%R;
}

sub checkInfo {
	my %r=%{$_[0]};
	my $rate=$_[1];
	my $minCov=$_[2];
	my @I=split(/\,/,$r{"I16"});
	my @QS=split(/\,/,$r{"QS"});
	if($r{"DP"}<$minCov){
		return 0;
	}
	my $s=$I[2]+$I[3];
#	if (($s/$r{"DP"})>$rate){
	if ($QS[1]>=0.001){
		return 1;
	}
	return 0;
}

sub parseInfo {
	my $info=shift;
	my @info=split(/\;/,$info);
	my %h;
	foreach my $stanza (@info){
		my ($key,$value)=split(/\=/,$stanza);
		$h{$key}=$value;
	}
	return \%h;
}

sub collectTheGarbage {
	my @files = @_;
	foreach my $file (@files){
		my $command="rm -rf $file";
		next unless -e $file;
		warn $command."\n";
		`$command`;
	}
	return 1;
}

sub prepFinal {
	my $finalDir = shift @_;
	my @files = @_;
	foreach my $file (@files){
		my $sPath=$file;
		my $oPath=$file;
		$oPath=~s/.+\///g;
		$oPath=$finalDir."/".$oPath;
		my $command = "mv $sPath $oPath";
		warn $command."\n";
		`$command`;
	}
	return 1;
}

