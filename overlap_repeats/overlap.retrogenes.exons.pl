#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $refblocks=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $gene=$s[0];
	my $idexon=$s[1];
	my $chr=$s[2];
	my $start=$s[3]+0;
	my $end=$s[4]+0;
	my $strand=$s[5];

	if($strand eq "."){
	    $strand="NA";
	}

	if(exists $refblocks->{$gene}){
	    push(@{$refblocks->{$gene}{"start"}},$start);
	    push(@{$refblocks->{$gene}{"end"}},$end);
	}
	else{
	    $refblocks->{$gene}={"chr"=>$chr,"strand"=>$strand,"start"=>[$start],"end"=>[$end]};
	}
	    
	$line=<$input>;
    }
    close($input);
}

##########################################################################

sub orderExonBlocks{
    my $refblocks=$_[0];
    my $refordered=$_[1];
    
    my %hashstart;
    
    foreach my $gene (keys %{$refblocks}){
	my $nbblocks=@{$refblocks->{$gene}{"start"}};
	my $chr=$refblocks->{$gene}{"chr"};
	my $strand=$refblocks->{$gene}{"strand"};

	for(my $i=0;$i<$nbblocks;$i++){
	    my $tb=${$refblocks->{$gene}{"start"}}[$i];
	    my $te=${$refblocks->{$gene}{"end"}}[$i];
	
	    if(exists $hashstart{$chr}){
		if(exists $hashstart{$chr}{$tb}){
		    if(exists $hashstart{$chr}{$tb}{$te}){
			if(exists $hashstart{$chr}{$tb}{$te}{$strand}){
			    push(@{$hashstart{$chr}{$tb}{$te}{$strand}},$gene);
			}
			else{
			    $hashstart{$chr}{$tb}{$te}{$strand}=[$gene];
			}
		    }
		    else{
			$hashstart{$chr}{$tb}{$te}={$strand=>[$gene]};
		    }
		}
		else{
		    $hashstart{$chr}{$tb}={$te=>{$strand=>[$gene]}};
		}
	    }
	    else{
		$hashstart{$chr}={$tb=>{$te=>{$strand=>[$gene]}}};
	    }
	}
	
    }

    foreach my $chr (keys %hashstart){
	$refordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"gene"=>[]};

	my @uniquestart=keys %{$hashstart{$chr}};
	
	my @sortedstart = sort {$a <=> $b} @uniquestart;

	foreach my $start (@sortedstart){
	    my @uniqueend=keys %{$hashstart{$chr}{$start}};
	    my @sortedend = sort {$a <=> $b} @uniqueend;

	    foreach my $end (@sortedend){
		foreach my $strand (keys %{$hashstart{$chr}{$start}{$end}}){
		    foreach my $gene (@{$hashstart{$chr}{$start}{$end}{$strand}}){
			push(@{$refordered->{$chr}{"start"}},$start);
			push(@{$refordered->{$chr}{"end"}},$end);
			push(@{$refordered->{$chr}{"strand"}},$strand);
			push(@{$refordered->{$chr}{"gene"}},$gene);
		    }
		} 
	    }
	} 
    }

}



##############################################################

sub readRetrogenes{
    my $pathin=$_[0];
    my $retrogenes=$_[1];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    my $line=<$input>;
   
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[0];
	my $prefix=substr $chr, 0,3;

	if($prefix eq "chr"){
	    $chr=substr $chr,3;
	}

	my $start=$s[1]+0;
	my $end=$s[2]+0;

	my $ss=$s[5];
	my $strand="NA";
	
	if($ss eq "+"){
	    $strand="1";
	} else{
	    if($ss eq "-"){
		$strand="-1";
	    } else{
		print "Weird strand!\n";
		print $line."\n";
		exit(1);
	    } 
	}
	
	my $id=$chr.",".$start.",".$end.",".$strand;

	$retrogenes->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	
	$line=<$input>;
    }

    close($input);

}

##############################################################

sub orderCoords{
    my $exons=$_[0];
    my $refordered=$_[1];

    my %hashstart;
    
    foreach my $exid (keys %{$exons}){
	my $chr=$exons->{$exid}{"chr"};
	my $b=$exons->{$exid}{"start"};
	my $e=$exons->{$exid}{"end"};
	my $s=$exons->{$exid}{"strand"};
	
	if(exists $hashstart{$chr}){
	    if(exists $hashstart{$chr}{$b}){
		push(@{$hashstart{$chr}{$b}{"end"}},$e);
		push(@{$hashstart{$chr}{$b}{"strand"}},$s);
		push(@{$hashstart{$chr}{$b}{"id"}},$exid);
	    }
	    else{
		$hashstart{$chr}{$b}={"end"=>[$e],"strand"=>[$s],"id"=>[$exid]};
	    }
	}
	else{
	    $hashstart{$chr}={$b=>{"end"=>[$e],"strand"=>[$s],"id"=>[$exid]}};
	}
    }
    
    foreach my $chr (keys %hashstart){
	$refordered->{$chr}={"start"=>[], "end"=>[], "id"=>[], "strand"=>[]};
	
	my @uniquestart=keys %{$hashstart{$chr}};
	
	my @sortedstart = sort {$a <=> $b} @uniquestart;
	
	foreach my $b (@sortedstart){
	    
	    my $nbblocks=@{$hashstart{$chr}{$b}{"end"}};
	    
	    for(my $i=0;$i<$nbblocks;$i++){
		my $strand=${$hashstart{$chr}{$b}{"strand"}}[$i];
		push(@{$refordered->{$chr}{"start"}},$b);
		push(@{$refordered->{$chr}{"strand"}},$strand);
		push(@{$refordered->{$chr}{"end"}},${$hashstart{$chr}{$b}{"end"}}[$i]);
		push(@{$refordered->{$chr}{"id"}},${$hashstart{$chr}{$b}{"id"}}[$i]);
	    }
	}
    }
}

##########################################################################

sub overlapBlocks{
    my $refblocks1=$_[0]; ## ordered exon blocks
    my $refblocks2=$_[1]; ## ordered exon blocks
    my $type=$_[2];
    my $refoverlap=$_[3]; ## overlap with blocks
   
    foreach my $chr (keys %{$refblocks1}){
	if(exists $refblocks2->{$chr}){
	    my $nbblocks1=@{$refblocks1->{$chr}{"start"}};
	    my $nbblocks2=@{$refblocks2->{$chr}{"start"}};
	    
	    my $firstindex2=0;  ## this is where we start looking for overlap
	    
	    for(my $i=0; $i<$nbblocks1; $i++){
		my $start1=${$refblocks1->{$chr}{"start"}}[$i];
		my $end1=${$refblocks1->{$chr}{"end"}}[$i];
		my $gene1=${$refblocks1->{$chr}{"gene"}}[$i];
		my $strand1=${$refblocks1->{$chr}{"strand"}}[$i];

		my $j=$firstindex2;

		while($j<$nbblocks2 && ${$refblocks2->{$chr}{"end"}}[$j]<$start1){ ## there cannnot be any overlap before that 
		    $j++;
		}

		$firstindex2=$j;

		while($j<$nbblocks2 && ${$refblocks2->{$chr}{"start"}}[$j]<=$end1){  ## we stop looking for overlap if the start coordinate of the second set of blocks is larger than the end coordinate of block
		 
		    my $M=max($start1,${$refblocks2->{$chr}{"start"}}[$j]);
		    my $m=min($end1,${$refblocks2->{$chr}{"end"}}[$j]);
		    		   
		    my $strand2=${$refblocks2->{$chr}{"strand"}}[$j];
		    
		    my $lenov=$m-$M+1;
		    
		    if($lenov>=1){

			if(($strand1 eq $strand2 && $type eq "sense") || ($strand1 ne $strand2 && $type eq "antisense") || ($type eq "both")){
			    my $keyblock=$start1.",".$end1;
			    
			    my $keyblock2=${$refblocks2->{$chr}{"start"}}[$j]." ".${$refblocks2->{$chr}{"end"}}[$j];
			    
			    my $frov1=($m-$M+1.0)/($end1-$start1+1.0);
			    my $frov2=($m-$M+1.0)/(${$refblocks2->{$chr}{"end"}}[$j]-${$refblocks2->{$chr}{"start"}}[$j]+1.0);
			    
			    if(exists $refoverlap->{$gene1}){
				if(!(exists $refoverlap->{$gene1}{$keyblock})){
				    $refoverlap->{$gene1}{$keyblock}={"start"=>[$M],"end"=>[$m]};
				}
				else{
				    push(@{$refoverlap->{$gene1}{$keyblock}{"start"}},$M);
				    push(@{$refoverlap->{$gene1}{$keyblock}{"end"}},$m);
				}
			    }
			    else{
				$refoverlap->{$gene1}={$keyblock=>{"start"=>[$M],"end"=>[$m]}};
			    }
			}
		    }
		    
		    
		    $j++;
		}
	    }
	}
    }
}

##########################################################################

sub makeOverlapBlocks{
    my $refoverlap=$_[0];
    my $refblocks=$_[1];

    foreach my $gene (keys %{$refoverlap}){

	$refblocks->{$gene}={};
	
	foreach my $idexon (keys %{$refoverlap->{$gene}}){

	    $refblocks->{$gene}{$idexon}={"start"=>[],"end"=>[]};

	    my %hashoverlap;
	    
	    my $nbpieces=@{$refoverlap->{$gene}{$idexon}{"start"}};
	    
	    for(my $i=0;$i<$nbpieces;$i++){
		my $b=${$refoverlap->{$gene}{$idexon}{"start"}}[$i];
		my $e=${$refoverlap->{$gene}{$idexon}{"end"}}[$i];

		if(exists $hashoverlap{$b}){
		    if($e>$hashoverlap{$b}){
			$hashoverlap{$b}=$e;
		    }
		}
		else{
		    $hashoverlap{$b}=$e;
		}
	    }

	    my @startoverlap;
	    my @endoverlap;

	    ## sort reference

	    my @uniquedb=keys %hashoverlap;
	    my @sorteddb= sort {$a <=> $b} @uniquedb;

	    ## update exon coordinates

	    my @startoverlap;
	    my @endoverlap;
	    	    
	    my $currentstart=$sorteddb[0];
	    my $currentend=$hashoverlap{$sorteddb[0]};
	    
	    for(my $k=1;$k<@sorteddb;$k++){
		my $thisstart=$sorteddb[$k];
		my $thisend=$hashoverlap{$sorteddb[$k]};
	    
		if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		    if($thisend>$currentend){
			$currentend=$thisend;
		    }
		}
		else{
		    push(@startoverlap,$currentstart);
		    push(@endoverlap,$currentend);
		    
		    # print "merging ".$currentstart ." to ".$currentend."\n";
		    
		    $currentstart=$thisstart;
		    $currentend=$thisend;
		}
	    }
	    
	    ## don't forget the last block

	    push(@startoverlap,$currentstart);
	    push(@endoverlap,$currentend);
	    
	    push(@{$refblocks->{$gene}{$idexon}{"start"}},@startoverlap);
	    push(@{$refblocks->{$gene}{$idexon}{"end"}},@endoverlap);
	}
    }
}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts overlaps between exons and retrogenes. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

my %parameters;

$parameters{"pathExonBlocks"}="NA";
$parameters{"pathRetrogenes"}="NA";
$parameters{"type"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonBlocks", "pathRetrogenes", "type","pathOutput");

my %numericpars;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## check if help was asked 

foreach my $arg (@ARGV){
    if($arg eq "--help"){
	printHelp(\@defaultpars, \%defaultvalues);
	exit(0);
    }
}

## check new parameters

my $nbargs=@ARGV;

for(my $i=0; $i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
	
	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0.0;
	}
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}

## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

##########################################################################
##########################################################################

print "Reading exon blocks...\n";
my %exons;
readExonBlocks($parameters{"pathExonBlocks"},\%exons);
print "Done\n";

print "Ordering exons...\n";
my %orderedexons;
orderExonBlocks(\%exons,\%orderedexons);
print "Done.\n";

##########################################################################

print "Reading retrogenes...\n";

my %retrogenes;
   
my @paths=split(",", $parameters{"pathRetrogenes"});

foreach my $path (@paths){
    print "Reading from ".$path."\n";
    readRetrogenes($path, \%retrogenes);
}

my %orderedretro;

orderCoords(\%retrogenes, \%orderedretro);

print "Done.\n";

##########################################################################

print "Extracting overlaps...\n";

my %overlap;
my $type=$parameters{"type"};

print "Overlap type: ".$type."\n";

overlapBlocks(\%orderedexons,\%orderedretro, $type, \%overlap); 

print "Done.\n";

##########################################################################

print "Making overlap blocks...\n";

my %ovblocks;

makeOverlapBlocks(\%overlap,  \%ovblocks);

print "Done.\n";

##########################################################################

print "Writing output...\n";

open(my $output,">".$parameters{"pathOutput"});

my $line="GeneID\tTotalLength\tOverlapLength";

print $output $line."\n";

foreach my $gene (keys %exons){
    my $chr=$exons{$gene}{"chr"};
    my $nbexons=@{$exons{$gene}{"start"}};
    
    my $totlength=0;

    for(my $i=0;$i<$nbexons;$i++){
	$totlength+=(${$exons{$gene}{"end"}}[$i]-${$exons{$gene}{"start"}}[$i]+1);
    }

    my $line=$gene."\t".$totlength;

    my $overlaplength=0;
    
    if(exists $ovblocks{$gene}){	    
	foreach my $idexon (keys %{$ovblocks{$gene}}){
	    my $nbpieces=@{$ovblocks{$gene}{$idexon}{"start"}};
	    
	    for(my $i=0;$i<$nbpieces;$i++){
		$overlaplength+=(${$ovblocks{$gene}{$idexon}{"end"}}[$i]-${$ovblocks{$gene}{$idexon}{"start"}}[$i]+1);
	    }
	}
    }
    
    $line.="\t".$overlaplength;
    
    print $output $line."\n";
    
}
close($output);


print "Done.\n";

##########################################################################
