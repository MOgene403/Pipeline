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

die "usage: perl $0 <Config file>\n\n" unless $#ARGV==0;
my $Config = Configuration->new($ARGV[0]);

warn "Finding Vectors...\n";
my $DataDir 	= $Config->get("DIRECTORIES","filtered_dir");
my $OutDir  	= $Config->get("DIRECTORIES","output_dir");
my $RefDir	= $Config->get("DIRECTORIES","reference");
my $TempDir	= $Config->get("DIRECTORIES","temp_dir");
my $vecDir 	= $Config->get("DIRECTORIES","vector_dir");
my $samtools 	= $Config->get("PATHS","samtools");
my $bwa		= $Config->get("PATHS","bwa");
my $snpRate	= $Config->get("OPTIONS","snpRate");
my $minCov	= $Config->get("OPTIONS","minCov");

my @Lines = $Config->getAll("CELL_LINE");

our %AFasta;
foreach my $grp (@Lines){
	my $vector = $Config->get("VECTORS",$grp);
	my $bwaRef = $OutDir."/".$vector.".ref.fasta";
	if(defined($AFasta{$vector})){
	}else{
		warn "Parsing $bwaRef\n";
		$AFasta{$vector}=Tools->LoadFasta($bwaRef);
	}
}


foreach my $grp (@Lines){
	my $bwaRef	= $OutDir."/".$Config->get("VECTORS",$grp).".ref.fasta";
	my $base 	= $Config->get("CELL_LINE",$grp);
	my $vector 	= $Config->get("VECTORS",$grp);
	my $outputDir 	= $OutDir."/".$base;
	my $bwaRoot	= $outputDir."/$base.Alignments";
	my $ins 	= $Config->get("INSERTS",$grp);

	opendir(VCF,$outputDir) || die "cannot open directory: $outputDir!\n$!\nexiting...\n";
	my @files=grep {m/vcf$/} readdir VCF;
	closedir VCF;
	my $j = $grp;
	my $insert = $Config->get("INSERTS",$grp);
	my %hash;
	foreach my $file (@files){
		my $path = $outputDir."/".$file;
		warn $path."\n";
		my @file = @{Tools->LoadFile($path)};
		foreach my $line (@file){
			next if $line=~m/^\#/;
			my @line=split(/\t/,$line);
			next unless(($line[0] eq $insert) || ($line[7]=~m/$insert/));
	#		print $line."\n";
			$hash{$line[7]}=$line;		
#21119   162     TRA00000003     N       <TRA>   .       PASS    IMPRECISE;CIEND=-708,708;CIPOS=-708,708;SVTYPE=TRA;SVMETHOD=EMBL.DELLYv0.5.4;CHR2=10;END=10386882;SVLEN=0;CT=5to5;PE=52;MAPQ=60 GT:GL:GQ:FT:RC:DR:DV:RR:RV       1/1:-308.7,-15.6534,0:157:PASS:0:0:52:0:0
#21119   15852   TRA00000001     N       <TRA>   .       LowQual IMPRECISE;CIEND=-489,489;CIPOS=-489,489;SVTYPE=TRA;SVMETHOD=EMBL.DELLYv0.5.4;CHR2=10;END=118244057;SVLEN=0;CT=3to3;PE=2;MAPQ=33 GT:GL:GQ:FT:RC:DR:DV:RR:RV       1/1:-195.996,-11.1346,0:111:PASS:0:0:37:0:0
#21119   15881   TRA00000002     N       <TRA>   .       PASS    IMPRECISE;CIEND=-263,263;CIPOS=-263,263;SVTYPE=TRA;SVMETHOD=EMBL.DELLYv0.5.4;CHR2=10;END=10386940;SVLEN=0;CT=3to3;PE=32;MAPQ=60 GT:GL:GQ:FT:RC:DR:DV:RR:RV       1/1:-195.996,-11.1346,0:111:PASS:0:0:37:0:0
		}
	}
	my @k = keys %hash;
	parseVCFtoCSV(\%hash,$j,$vector);	
}


sub parseVCFtoCSV {
	my %R = %{$_[0]};
	my $j=$_[1];
	my $vector = $_[2];
	print "Source Sample,Contig A,A_coordinate,Type,Contig B,B_coordinate,Length,Orientation,NumReads,Segment A,Segment B\n";
	foreach my $key (keys %R){
		my $line=$R{$key};
		my @line=split(/\t/,$line);
		my %info = %{parseInfo($line[7],$line)};
		$line[2] =~ s/\d+//;
		my $nLine = $j.",".$line[0].",".$line[1].",".$line[2];
		$nLine .= ",".$info{CHR2};
		$nLine .= ",".$info{END};
		$nLine .= ",".$info{SVLEN};
		$nLine .= ",".$info{CT};
		$nLine .= ",".$info{PE};
		my ($segA,$segB)=split(/\,/,getSegment($line[2],$line[0],$line[1],$info{END},$info{CHR2},$info{CT},$vector));
		$nLine .=",".$segA;
		$nLine .=",".$segB;
		print $nLine."\n";
	}
}

sub getSegment { 
	my $type=shift;
	my $chrA = shift;
	my $coordA = shift;
	my $coordB = shift;
	my $chrB = shift;
	my $ct = shift;
	my $vector = shift;
	my %Fasta=%{$AFasta{$vector}};
	my $segmentA="";
	my $segmentB="";
	if($type=~m/TRA/i){
		if(defined($Fasta{$chrA})){
			my $start = $coordA - 100;
			$start =0  if $start < 0;
			my $length = 200;
			$length = length($Fasta{$chrA}) - $coordA if ($length+$coordA)>length($Fasta{$chrA});
			$segmentA = substr($Fasta{$chrA},$start,$length);
		}else{
			die "cannot find chromosome $chrA!\n";
		}
		if(defined($Fasta{$chrB})){
			my $start = $coordB - 100;
			$start =0  if $start < 0;
			my $length = 200;
			$length = length($Fasta{$chrB}) - $coordB if ($length+$coordB)>length($Fasta{$chrB});
			$segmentB = substr($Fasta{$chrB},$start,$length);
		}else{
			die "cannot find chromosome $chrA!\n";
		}
	}elsif($type=~m/DEL/i){
		if(defined($Fasta{$chrA})){
			my $start = $coordA - 100;
			$start =0  if $start < 0;
			my $length = 200;
			$length = length($Fasta{$chrA}) - $coordA if ($length+$coordA)>length($Fasta{$chrA});
			$segmentA = substr($Fasta{$chrA},$start,$length);
		}else{
			die "cannot find chromosome $chrA!\n";
		}
		if((abs($coordB-$coordA))>200){
			my $start = $coordB - 100;
			$start =0  if $start < 0;
			my $length = 200;
			$length = length($Fasta{$chrA}) - $coordA if ($length+$coordA)>length($Fasta{$chrA});
			$segmentB = substr($Fasta{$chrA},$start,$length);
			
		}
	}else{
	}
	return $segmentA.",".$segmentB;

}


sub parseInfo {
	my $info=shift;
	my $line=shift;
	my %H;
	my @I = split(/\;/,$info);
	foreach my $i (@I){
		my ($k,$v)=split(/\=/,$i);
		$H{$k}=$v;
	}
#SNP VCF
#21119   15862   .       A       X       0       .       DP=70;I16=0,70,0,0,2294,77842,0,0,3405,168105,0,0,237,1779,0,0;QS=1.000000,0.000000,0.000000,0.000000       PL      0,211,216
#21119   15863   .       A       X       0       .       DP=42;I16=0,42,0,0,1504,54694,0,0,2024,99644,0,0,195,1347,0,0;QS=1.000000,0.000000,0.000000,0.000000        PL      0,126,206
#21119   15864   .       T       X       0       .       DP=42;I16=0,42,0,0,1412,49384,0,0,2024,99644,0,0,153,999,0,0;QS=1.000000,0.000000,0.000000,0.000000 PL      0,126,203
#21119   15865   .       A       X       0       .       DP=39;I16=0,39,0,0,1339,47855,0,0,1874,92144,0,0,114,732,0,0;QS=1.000000,0.000000,0.000000,0.000000 PL      0,117,203
#21119   15866   .       T       X       0       .       DP=39;I16=0,39,0,0,1343,47799,0,0,1874,92144,0,0,75,543,0,0;QS=1.000000,0.000000,0.000000,0.000000  PL      0,117,203
#
#
#TRA VCF (delly) :
#21119   132     TRA00000002     N       <TRA>   .       PASS    IMPRECISE;CIEND=-487,487;CIPOS=-487,487;SVTYPE=TRA;SVMETHOD=EMBL.DELLYv0.5.4;CHR2=3;END=18639948;SVLEN=0;CT=5to5;PE=93;MAPQ=60      GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-550.988,-30.9945,0:310:PASS:0:0:103:0:0
#21119   132     TRA00000001     N       <TRA>   .       PASS    IMPRECISE;CIEND=-584,584;CIPOS=-584,584;SVTYPE=TRA;SVMETHOD=EMBL.DELLYv0.5.4;CHR2=3;END=39196983;SVLEN=0;CT=5to3;PE=6;MAPQ=43       GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-550.988,-30.9945,0:310:PASS:0:0:103:0:0
#21119   15835   TRA00000003     N       <TRA>   .       LowQual IMPRECISE;CIEND=-486,486;CIPOS=-486,486;SVTYPE=TRA;SVMETHOD=EMBL.DELLYv0.5.4;CHR2=3;END=18639996;SVLEN=0;CT=3to3;PE=2;MAPQ=45       GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-10.5972,-0.900273,0:9:LowQual:0:0:3:0:0
#
	$H{CHR2}="NA" unless defined $H{CHR2};
	$H{END}="NA" unless defined $H{END};
	$H{SVLEN}="0" unless defined $H{SVLEN};
	$H{CT}="0" unless defined $H{CT};
	if(defined($H{DP})){
		$H{PE} = $H{DP};
	}
	$H{PE}="0" unless defined $H{PE};
	#check(\%H,$line);
	return \%H;
}

sub check {
	my %hash=%{$_[0]};
	my $line=$_[1];
	die "malformed vcf\n$line\nNo:CHR2 \n" unless(defined($hash{CHR2}));
	die "malformed vcf\n$line\nNo:SVLEN \n" unless(defined($hash{SVLEN}));
	die "malformed vcf\n$line\nNo:END \n" unless(defined($hash{END}));
	die "malformed vcf\n$line\nNo:CT \n" unless(defined($hash{CT}));
	die "malformed vcf\n$line\nNo:PE \n" unless(defined($hash{PE}));
	die "malformed vcf\n$line\nNo:MAPQ \n" unless(defined($hash{MAPQ}));
	die "malformed vcf\n$line\nNo:CIEND \n" unless(defined($hash{CIEND}));
	die "malformed vcf\n$line\nNo:CIPOS \n" unless(defined($hash{CIPOS}));
	return 1;
}

#[INSERTS]
#FSR14004-1508_S6 = 17629
#FSR14004-1515_S7 = 20791
#FSR14004-B10-123_S34 = 21119
#FSR14004-B12-124_S36 = 21119
#FSR14004-C02-131_S38 = 21119
#FSR14004-C03-126_S39 = 21119


exit(0);


