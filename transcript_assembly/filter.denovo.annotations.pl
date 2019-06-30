#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

########################################################################

sub readReadthroughTranscripts{
    my $pathin=$_[0];
    my $rt=$_[1];
    
    open(my $input, $pathin);
    my $line=<$input>; ## header
    $line=<$input>; ## first actual line
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $txid=$s[0];

	$rt->{$txid}=1;

	$line=<$input>;
    }

    close($input);
}

########################################################################

sub readSpliceJunctions{
    my $pathin=$_[0];
    my $junctions=$_[1];
    
    open(my $input, $pathin);
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>; ## first actual line
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $txid=$s[$header{"TranscriptID"}];
	my $nbintrons=$s[$header{"NbIntrons"}]+0;
	my $nbok=$s[$header{"NbSupportedIntrons"}]+0;
	my $nbwrong=$s[$header{"NbWrongStrand"}]+0;
	
	$junctions->{$txid}={"nbintrons"=>$nbintrons, "nbok"=>$nbok,"nbwrong"=>$nbwrong};

	$line=<$input>;
    }

    close($input);

}

########################################################################

sub readEnsemblSynonyms{
    my $pathin=$_[0];
    my $synonyms=$_[1];
    
    open(my $input, $pathin);
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>; ## first actual line
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $txid=$s[$header{"DeNovoID"}];
	my $ensid=$s[$header{"EnsemblID"}];
	
	if(exists $synonyms->{$txid}){
	    print "Weird! ".$txid." appears more than once in the synonyms file.\n";
	    $synonyms->{$txid}{$ensid}=1;
	}
	else{
	   $synonyms->{$txid}={$ensid=>1};
	}
	
	$line=<$input>;
    }

    close($input);
}

########################################################################

sub readCoverage{
    my $pathin=$_[0];
    my $coverage=$_[1];
    
    open(my $input, $pathin);
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>; ## first actual line
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $txid=$s[$header{"TranscriptID"}];
	my $covsense=$s[$header{"CoverageSense"}]+0.0;
	my $covantisense=$s[$header{"CoverageAntisense"}]+0.0;
	
	$coverage->{$txid}={"sense"=>$covsense, "antisense"=>$covantisense};

	$line=<$input>;
    }

    close($input);
}

########################################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];
    
    my $res="NA";
    
    my @grepres=grep(/${pattern}/,@{$array});

    my $nbg=@grepres;

    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    }
    
    return $res;
}

###########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script filters de novo annotations.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################
##########################################################################

## parameters 

my %parameters;
$parameters{"pathAssembledGTF"}="NA";
$parameters{"pathEnsemblSynonyms"}="NA";
$parameters{"pathReadthroughTranscripts"}="NA";
$parameters{"pathCoverage"}="NA";
$parameters{"minRatioSenseAntisense"}="NA";
$parameters{"pathJunctions"}="NA";
$parameters{"pathOutputGTF"}="NA";
$parameters{"pathOutputDiscardedTranscripts"}="NA";

my @defaultpars=("pathAssembledGTF", "pathEnsemblSynonyms", "pathReadthroughTranscripts", "pathCoverage", "minRatioSenseAntisense", "pathJunctions", "pathOutputGTF", "pathOutputDiscardedTranscripts");

my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
	
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
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

##############################################################
##############################################################

print "Reading readthrough transcript ids...\n";
my %readthrough;

readReadthroughTranscripts($parameters{"pathReadthroughTranscripts"}, \%readthrough);

my $nbtx=keys %readthrough;

print "Found ".$nbtx." read-through transcripts.\n";

print "Done.\n";

##############################################################

print "Reading splice junctions info...\n";

my %junctions;
readSpliceJunctions($parameters{"pathJunctions"}, \%junctions);

my $nbtx=keys %junctions;

print "Found ".$nbtx." transcripts with splice junction information.\n";

print "Done.\n";

##############################################################

print "Reading sense/antisense transcript coverage...\n";

my %coverage;
readCoverage($parameters{"pathCoverage"}, \%coverage);

my $nbtx=keys %coverage;

print "Found ".$nbtx." transcripts with coverage information.\n";

print "Done.\n";

##############################################################

print "Reading Ensembl synonyms...\n";

my %synonyms;
readEnsemblSynonyms($parameters{"pathEnsemblSynonyms"}, \%synonyms);

my $nbtx=keys %synonyms;

print "Found ".$nbtx." transcripts with Ensembl synonyms.\n";

print "Done.\n";

##############################################################

print "Reading and filtering GTF...\n";

my %discarded;
my %kept;

open(my $outputgtf, ">".$parameters{"pathOutputGTF"});
open(my $input, $parameters{"pathAssembledGTF"});

my $minratiosas=$parameters{"minRatioSenseAntisense"}+0.0;

print "We keep transcripts if their ratio sense/antisense is at least ".$minratiosas.".\n";

my $line=<$input>;
my %unstranded;

while($line){
    chomp $line;
    my @s=split("\t", $line);
	
    my $type=$s[2];

    if($type eq "exon"){
	## we only write exon lines
	
	my $info=$s[8];
	my @infoarray=split(";", $info);
	
	my $txid=findInfo("transcript_id", \@infoarray);
	my $prefix=substr $txid,0,3;

	
	my $ss=$s[6];
	my $strand="NA";
	
	if($ss eq "+"){
	    $strand="1";
	} else{
	    if($ss eq "-"){
		$strand="-1";
	    } 
	}
	
	if($strand ne "NA"){
	    if(exists $synonyms{$txid}){
		$kept{$txid}="Ensembl_synonym";
		print $outputgtf $line."\n";
	    }
	    else{
		if($prefix eq "ENS"){
		    $kept{$txid}="Ensembl";
		    print $outputgtf $line."\n";
		}
		else{
		    ## purely de novo
		    
		    if((!exists $coverage{$txid}) || (!exists $junctions{$txid})){
			print "Weird!!! cannot find ".$txid." in coverage or junctions input.\n";
			exit(1);
		    }
		    
		    if(exists $readthrough{$txid}){
			$discarded{$txid}="readthrough";
		    }
		    else{
			my $covsense=$coverage{$txid}{"sense"};
			my $covanti=$coverage{$txid}{"antisense"};
			
			if($covsense==0){
			    $discarded{$txid}="nocoverage";
			}
			else{
			    my $thisratio=1;
			    
			    if($covanti>0){
				$thisratio=($covsense+0.0)/($covanti+0.0);
			    }
			    
			    if($thisratio < $minratiosas){
				$discarded{$txid}="lowSAS";
			    }
			    else{
				if($junctions{$txid}{"nbwrong"}>0){
				    $discarded{$txid}="wrongstrand";
				}
				else{
				    ## everything ok
				    print $outputgtf $line."\n";
				}
			    }
			}
		    }
		}
	    }
	}
	else{
	    $discarded{$txid}="unstranded";
	    $unstranded{$txid}=1;
	}
    }

    $line=<$input>;
}

close($outputgtf);
close($input);

my $nbunstranded=keys %unstranded;

print "Removed ".$nbunstranded." unstranded transcripts.\n";

print "Done.\n";

##############################################################

print "Writing output for discarded transcripts...\n";

my %reasons;

foreach my $txid (keys %discarded){
    my $r=$discarded{$txid};
    
    if(exists $reasons{$r}){
	$reasons{$r}++;
    }
    else{
	$reasons{$r}=1;
    }
}

foreach my $r (keys %reasons){
    print $reasons{$r}." discarded: ".$r."\n";
}

##############################################################

open(my $outputdiscarded, ">".$parameters{"pathOutputDiscardedTranscripts"});
print $outputdiscarded "TranscriptID\tReason\n";
 
foreach my $txid (keys %discarded){
    print $outputdiscarded $txid."\t".$discarded{$txid}."\n";
}

close($outputdiscarded);

print "Done.\n";

##############################################################
