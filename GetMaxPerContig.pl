#!/usr/bin/perl
use warnings;
use strict;

my %Maxes;
while(<>){
	chomp $_;
	my ($contig,$pos,$depth)=split(/\t/,$_);
	if(defined($Maxes{$contig})){
		if($depth>$Maxes{$contig}){
			$Maxes{$contig}=$depth;
		}
	}else{
		$Maxes{$contig}=$depth;
	}
}

map {print $_."\t".$Maxes{$_}."\n"} sort {$Maxes{$b} <=> $Maxes{$a}} keys %Maxes;
