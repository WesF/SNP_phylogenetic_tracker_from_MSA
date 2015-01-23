#!/usr/bin/perl -w

use warnings;

my $INHANDLE = $ARGV[0];
open(IN,"$INHANDLE")|| die "cannot open $INHANDLE";


my @Other = (0,0,0,0,0,0,0,0,0,0,0,0);
my %phylastats=(Firmicutes=>[0,0,0,0,0,0,0,0,0,0,0,0],
				'Bacteroidetes/Chlorobi group'=>[0,0,0,0,0,0,0,0,0,0,0,0],
				Actinobacteria=>[0,0,0,0,0,0,0,0,0,0,0,0],
				Proteobacteria=>[0,0,0,0,0,0,0,0,0,0,0,0],
				Fusobacteria=>[0,0,0,0,0,0,0,0,0,0,0,0]
				);
				
while(<IN>){
	my @terms = split(/\t/, $_);
	if(exists $phylastats{$terms[12]}){
		if($terms[1] eq 'N' ){
			$phylastats{$terms[12]}[0]++;
		}
		if($terms[2] eq 'P' ){
			$phylastats{$terms[12]}[1]++;
		}
		if($terms[3] eq 'F' or $terms[3] eq 'Y' or $terms[3] eq 'P'){
			$phylastats{$terms[12]}[2]++;
		}
		if($terms[4] eq 'L' or $terms[4] eq 'H'){
			$phylastats{$terms[12]}[3]++;
		}
		if($terms[5] eq 'G' or $terms[5] eq 'Y'){
			$phylastats{$terms[12]}[4]++;
		}
		if($terms[6] eq 'F' or $terms[6] eq 'Y'){
			$phylastats{$terms[12]}[5]++;
		}
		if($terms[7] eq 'Y' or $terms[7] eq 'L' or $terms[7] eq 'P'){
			$phylastats{$terms[12]}[6]++;
		}
		if($terms[8] eq 'S' or $terms[8] eq 'L' or $terms[8] eq 'H'){
			$phylastats{$terms[12]}[7]++;
		}
		if($terms[9] eq 'F' ){
			$phylastats{$terms[12]}[8]++;
		}
		if($terms[10] eq 'L'){
			$phylastats{$terms[12]}[9]++;
		}
		if($terms[11] eq 'F'){
			$phylastats{$terms[12]}[10]++;
		}
		if($terms[3] ne '-' ){
			$phylastats{$terms[12]}[11]++;
		}
	}else{
	
		
		if($terms[1] eq 'N' ){
			$Other[0]++;
		}
		if($terms[2] eq 'P' ){
			$Other[1]++;
		}
		if($terms[3] eq 'F' or $terms[3] eq 'Y' or $terms[3] eq 'P'){
			$Other[2]++;
		}
		if($terms[4] eq 'L' or $terms[4] eq 'H'){
			$Other[3]++;
		}
		if($terms[5] eq 'G' or $terms[5] eq 'Y'){
			$Other[4]++;
		}
		if($terms[6] eq 'F' or $terms[6] eq 'Y'){
			$Other[5]++;
		}
		if($terms[7] eq 'Y' or $terms[7] eq 'L' or $terms[7] eq 'P'){
			$Other[6]++;
		}
		if($terms[8] eq 'S' or $terms[8] eq 'L' or $terms[8] eq 'H'){
			$Other[7]++;
		}
		if($terms[9] eq 'F' ){
			$Other[8]++;
		}
		if($terms[10] eq 'L'){
			$Other[9]++;
		}
		if($terms[11] eq 'F'){
			$Other[10]++;
		}
		if($terms[3] ne '-' ){
			$Other[11]++;
		}
		
		
	}
	
	
}

foreach my $key (keys %phylastats){
	print "$key\t$phylastats{$key}[0]\t$phylastats{$key}[1]\t$phylastats{$key}[2]\t$phylastats{$key}[3]\t$phylastats{$key}[4]\t$phylastats{$key}[5]\t$phylastats{$key}[6]\t$phylastats{$key}[7]\t$phylastats{$key}[8]\t$phylastats{$key}[9]\t$phylastats{$key}[10]\t$phylastats{$key}[11]\n";
	
}
print "Other\t$Other[0]\t$Other[1]\t$Other[2]\t$Other[3]\t$Other[4]\t$Other[5]\t$Other[6]\t$Other[7]\t$Other[8]\t$Other[9]\t$Other[10]\t$Other[11]\n";
close(IN);