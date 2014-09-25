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

my $ref=$Config->get("PATHS","reference_dir");
warn "Finding Vectors...\n";
my $vecDir = $Config->get("DIRECTORIES","vector_dir");
my @Lines = $Config->getAll("INSERTS");
my $VCFdir = $Config->get("DIRECTORIES","output_dir")."/VCFs";

opendir(VCF,$VCFdir) || die "cannot open directory: $VCFdir!\n$!\nexiting...\n";
my @files=grep {m/vcf$/} readdir VCF;
closedir VCF;

our %Fasta;
foreach my $line (@Lines){
	warn "parsing $line...\n";
	my $ref=$Config->get("DIRECTORIES","output_dir")."/".$Config->get("VECTORS",$line).".ref.fasta";
	%Fasta = %{Tools->LoadFasta($ref)};
	my $j = $line;
	my $insert = $Config->get("INSERTS",$line);
	my @files = grep {m/$line/} @files;
	print "\n";
#	print join("\n",@files)."\n";
	my %hash;
	foreach my $file (@files){
		my $path = $VCFdir."/".$file;
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
	parseVCFtoCSV(\%hash,$j);	
}


sub parseVCFtoCSV {
	my %R = %{$_[0]};
	my $j=$_[1];
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
		my ($segA,$segB)=split(/\,/,getSegment($line[2],$line[0],$line[1],$info{END},$info{CHR2},$info{CT}));
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
	check(\%H,$line);
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


