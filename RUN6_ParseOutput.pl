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

our @OUTPUT;
our @FS;
my $outputDir	= $OutDir."/VCFs";	
my $lbfs = $outputDir."/FS.fasta.txt";
my $output = $outputDir."/SNPs.csv";

foreach my $grp (@Lines){
	my $bwaRef	= $OutDir."/".$Config->get("VECTORS",$grp).".ref.fasta";
	my $exclusions	=$Config->get("DIRECTORIES","vector_dir")."/".$Config->get("VECTORS",$grp).".exclusions.tab";
	my %E = %{_LoadExclusions($exclusions,200,1000)};
	my $base 	= $Config->get("CELL_LINE",$grp);
	my $nk 		= $base;
	$nk 		=~ s/\_S\d+//;
	my $out 	= $Config->get("OUTPUT",$nk);
	my $OFile = "does not matter\n";
	my $vector 	= $Config->get("VECTORS",$grp);
	my $inputDir	= $OutDir."/".$base;
	my $bwaRoot	= $outputDir."/$base.Alignments";
	my $ins 	= $Config->get("INSERTS",$grp);
	opendir(VCF,$inputDir) || die "cannot open directory: $outputDir!\n$!\nexiting...\n";
	my @files=grep {m/vcf$/} readdir VCF;
	closedir VCF;
	my $j = $grp;
	my $insert = $Config->get("INSERTS",$grp);
	my %hash;
	foreach my $file (@files){
		my $path = $inputDir."/".$file;
		warn $path."\n";
		my @file = @{Tools->LoadFile($path)};
		foreach my $line (@file){
			next if $line=~m/^\#/;
			my @line=split(/\t/,$line);
			next unless(($line[0] eq $insert) || ($line[7]=~m/$insert/));
			$hash{$line[7]}=$line;		
		}
	}
	my @k = keys %hash;
	parseVCFtoCSV(\%hash,$j,$vector,$OFile,$out,\%E);
}


Tools->printToFile($output,\@OUTPUT);
Tools->printToFile($lbfs,\@FS);

sub _LoadExclusions {
	my $file=shift;
	my $minSize=shift;
	my $flank=shift;
	my @file=@{Tools->LoadFile($file)};
	my %E;
	foreach my $line (@file){
		chomp $line;
		my ($vect,$ref,$ident,$length,$mm,$gap,$qstart,$qend,$rstart,$rend,$eval,$bit)=split(/\t/,$line);
		next unless $length > $minSize;
		my $B=$rstart-$flank;
		my $E=$rend+$flank;
		warn "Loading Exclusion: $ref from $B to $E\n";
		if(defined($E{$ref})){
			for(my$i=$B;$i<=$E;$i++){
				$E{$ref}{$i}=1;
			}
		}else{
			$E{$ref}={};
			for(my$i=$B;$i<=$E;$i++){
				$E{$ref}{$i}=1;
			}
		}
	}
	return \%E;
#21119   9       99.73   2623    6       1       3900    6522    11500049        11502670        0.0      4804
#21119   9       99.17   1202    1       4       8497    9698    11506510        11507702        0.0      2156
#21119   4       80.61   196     20      9       13504   13687   163046119       163046308       5e-28     135
#21119   3       80.61   196     20      9       13504   13687   160830037       160829848       5e-28     135

}

sub parseVCFtoCSV {
	my %R = %{$_[0]};
	my $j=$_[1];
	my $vector = $_[2];
	my $File=$_[3];
	my $base = $_[4];
	my %Exc  = %{$_[5]};
	my @O;
	push @O, "Source Sample,Contig A,A_coordinate,Type,Contig B,B_coordinate,Orientation,NumReads,Sequence,(SNP) Ref Base,(SNP) Alt Base, (SNP) Ref Ratio, (SNP) Alt Ratio";
	my $rs=1;
	my $ls=1;
	foreach my $key (keys %R){
		my $line=$R{$key};
		my @line=split(/\t/,$line);
		my %info = %{parseInfo($line[7],$line)};
		$line[2] =~ s/\d+//;
		my $nLine = $base.",".$line[0].",".$line[1];
warn "WARNING HARDCODED FILTERING!!!!\n";
		if(defined($info{I16})){
		next unless ((($line[1] > 726)&&($line[1] < 3095)) || (($line[1] > 6529)&&($line[1] < 8475)) || (($line[1] > 11525)&&($line[1]<13468)) || (($line[1] > 14918)&&($line[1] < 15469)));
			if(my $cons=_checkI16($line[3],$line[4],\%info,$snpRate,$minCov)){
				$nLine .= ",SNP";
				$nLine .= ",NA";
				$nLine .= ",NA";
				$nLine .= ",NA";
				$nLine .= ",".$info{DP};
				$nLine .= ",$cons";
				push @OUTPUT, $nLine;
			}
		}else{
			next if defined $Exc{$info{CHR2}}{$info{END}};
			$nLine .= ",".$line[2];
			$nLine .= ",".$info{CHR2};
			$nLine .= ",".$info{END};
			$nLine .= ",".$info{CT};
			$nLine .= ",".$info{PE};
			my $sequence=getSegment($line[2],$line[0],$line[1],$info{END},$info{CHR2},$info{CT},$vector,200);
#			$nLine .=",".$sequence;
			if(($info{CT} eq "5to3")||($info{CT} eq "3to3")){ ## left
				push @FS, ">$base"."_LBFS-$ls\n$sequence";
				$ls++;	
			}else{ ## right
				push @FS, ">$base"."_RBFS-$rs\n$sequence";
				$rs++;	
			}			
			push @OUTPUT, $nLine;
		}
	}
#	push @O, "\n";
	#Tools->printToFile($File,\@O);
}

sub _checkI16 {
	my $cons = shift;
	my $Alts = shift;
	my %Info = %{$_[0]};
	my $rate = $_[1];
	my $min = $_[2];
	my $I16 = $Info{I16};
	my $QS  = $Info{QS};
	my @Alt = split(/\,/,$Alts);
	unshift @Alt, $cons;
	my @QS  = split(/\,/,$QS);
	my @I16=split(/\,/,$I16);
	if(($I16[0] == 0) && ($I16[1] == 0) && ($I16[2] == 0) && ($I16[3] ==0 )){
		return undef;
	}
	my $ret = "";
	my $freq = "";
	my $p=0;
	for(my $i=0;$i<=$#Alt;$i++){
		if($QS[$i] > $rate) {
			$ret.=",".$Alt[$i];
			$freq.= ",".$QS[$i];
			$p+=1;
		}
	}
	return $ret.",".$freq if $p > 1;
	return undef;
}

sub getSegment { 
	my $type=shift;
	my $chrA = shift;
	my $coordA = shift;
	my $coordB = shift;
	my $chrB = shift;
	my $ct = shift;
	my $vector = shift;
	my $N = shift;
	my %Fasta=%{$AFasta{$vector}};
	my $segmentA="";
	my $segmentB="";
	my $seq="";
	if($type=~m/TRA/i){
		if($ct eq "5to5"){
			## get first N coordinates of Chr A upstream of at Pos A, and first N coordinates of Chr B starting after Pos B
			my $start = $coordA - $N;
			$start = 0 if $start < 0;
			my $length = $N;
			$length = length($Fasta{$chrA}) - $coordA if ($length+$coordA)>length($Fasta{$chrA});
			$segmentA = substr($Fasta{$chrA},$start,$length);
			$segmentA = Tools->revcomp($segmentA);
			$start = $coordB;
			$length = $N;
			$length = length($Fasta{$chrB}) - $coordB if ($length+$coordB)>length($Fasta{$chrB});
			$segmentB = substr($Fasta{$chrB},$start,$length);
			$segmentA = lc($segmentA);
			$segmentB = uc($segmentB);
			$seq = $segmentA.$segmentB;
		}elsif($ct eq "5to3"){
			## Get first N nucleotides after A in ChrA, get first N nucleotides before B in chr B
			my $start = $coordA;
			my $length = $N;
			$length = length($Fasta{$chrA}) - $coordA if ($length+$coordA)>length($Fasta{$chrA});
			$segmentA = substr($Fasta{$chrA},$start,$length);
			$start = $coordB-$N;
			$start = 0 if $start < 0;
			$length = $N;
			$length = length($Fasta{$chrB}) - $coordB if ($length+$coordB)>length($Fasta{$chrB});
			$segmentB = substr($Fasta{$chrB},$start,$length);
			$segmentA = lc($segmentA);
			$segmentB = uc($segmentB);
			$seq = $segmentB.$segmentA;
			$seq = Tools->revcomp($seq);
		}elsif($ct eq "3to5"){
			## get first N nucleotides before A in ChrA, get first N nucleotides after B in Chr B
			my $start = $coordA-$N;
			$start = 0 if $start < 0;
			my $length = $N;
			$length = length($Fasta{$chrA}) - $coordA if ($length+$coordA)>length($Fasta{$chrA});
			$segmentA = substr($Fasta{$chrA},$start,$length);
			$start = $coordB;
			$length = $N;
			$length = length($Fasta{$chrB}) - $coordB if ($length+$coordB)>length($Fasta{$chrB});
			$segmentB = substr($Fasta{$chrB},$start,$length);
			$segmentA = lc($segmentA);
			$segmentB = uc($segmentB);
			$seq=$segmentA.$segmentB;
		}elsif($ct eq "3to3"){
			## get revcomp of the first N nucleotides upstream of A in Chr A, get N nucleotides upstream of B in Chr B
			my $start = $coordA - $N;
			$start = 0 if $start < 0;
			my $length = $N;
			$length = length($Fasta{$chrA}) - $coordA if ($length+$coordA)>length($Fasta{$chrA});
			$segmentA = substr($Fasta{$chrA},$start,$length);
			$segmentA = Tools->revcomp($segmentA);
			$start = $coordB-$N;
			$start = 0 if $start < 0;
			$length = $N;
			$length = length($Fasta{$chrB}) - $coordB if ($length+$coordB)>length($Fasta{$chrB});
			$segmentB = substr($Fasta{$chrB},$start,$length);
			$segmentA = lc($segmentA);
			$segmentB = uc($segmentB);
			$seq=$segmentB.$segmentA;
			$seq= Tools->revcomp($seq);
		}else{
		}
	}elsif($type=~m/DEL/i){
	}else{
	}
	return $seq;

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


