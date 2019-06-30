#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $exontx=$_[2];
    my $genetx=$_[3];
    my $txgene=$_[4];
    my $txex=$_[5];
  
    open(my $input, $pathin);
    
    my $line=<$input>;
    my $prefix=substr $line, 0,1;
    
    while($prefix eq "#"){
	$line=<$input>;
	$prefix=substr $line, 0,1;
    }

    my $nbunstranded=0; ## we discard genes with "." strand

    my %rttx;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon"){
	
	    my $chr=$s[0];
	    my $start=$s[3]+0; ## 1-based
	    my $end=$s[4]+0;
	    my $strand=$s[6];
	  

	    if($strand ne "."){
		my $info=$s[8];
		my @t=split(";", $info);
		my $geneid=findInfo("gene_id", \@t);
		my $txid=findInfo("transcript_id", \@t);
		
		if($txid eq "NA"){
		    print "could not find transcript in ".$line."\n";
		    exit(1);
		}		

		if($geneid eq "NA"){
		    print "could not find gene in ".$line."\n";
		    exit(1);
		}		
		
		my $exonid=$chr.",".$start.",".$end.",".$strand;
		
		## fill in exon coords
		
		$exoncoords->{$exonid}={"chr"=>$chr,"start"=>$start, "end"=>$end, "strand"=>$strand};
		
		## exon - tx correspondence
		
		if(exists $exontx->{$exonid}){
		    $exontx->{$exonid}{$txid}=1;
		}
		else{
		    $exontx->{$exonid}={$txid=>1};
		}
		
		if(exists $genetx->{$geneid}){
		    $genetx->{$geneid}{$txid}=1;
		}
		else{
		    $genetx->{$geneid}={$txid=>1};
		}
		
		if(exists $txgene->{$txid}){
		    if($txgene->{$txid} ne $geneid){
			print "Weird! ".$txid." is associated to more than one gene: ".$geneid. " ".$txgene->{$txid}."\n";
		    }
		} else{
		    $txgene->{$txid}=$geneid;
		}
		
		if(exists $txex->{$txid}){
		    $txex->{$txid}{$exonid}=1;
		}
		else{
		    $txex->{$txid}={$exonid=>1};
		}
	    }
	    else{
		$nbunstranded++;
	    }
	}
   	
	$line=<$input>;
    }
    
    close($input);

    print "We discarded ".$nbunstranded." exons with undefined strand.\n";
  
}

##############################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];
    
    my $res="NA";
    
    my @grepres=grep(/${pattern}/,@{$array});

    my $nbg=@grepres;
    my $nbreal=0;

    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    } else{
	my $nbreal=0;
	
	foreach my $g (@grepres){
	    $g =~ s/^\s+//; ## remove whitespace
	    my @u=split(" ",$g);

	    if($u[0] eq $pattern){
		$nbreal++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }
    
    if($nbreal>1){
	return "NA";
    }
    return $res;
}

##########################################################################

sub readSelectedTranscripts{
    my $pathin=$_[0];
    my $seltx=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $txid=$s[$header{"TranscriptID"}];
	my $geneid=$s[$header{"GeneID"}];
	my $decision=$s[$header{"Decision"}];
	my $reason=$s[$header{"Reason"}];
	my $ensid=$s[$header{"OverlappingEnsemblGenes"}];

	if($decision eq "Accepted"){
	    if($ensid eq "NA"){
		if($reason ne "NewGene"){
		    print "No Ensembl ID, but reason is not new gene.\n";
		    print $line."\n";
		    exit(1);
		}

		if(exists $seltx->{$geneid}){
		    $seltx->{$geneid}{$txid}=1;
		} else{
		    $seltx->{$geneid}={$txid=>1};
		}
	    } else{
		if($reason eq "NewGene"){
		    print "Existing Ensembl ID, but reason new gene.\n";
		    print $line."\n";
		    exit(1);
		}
		
		if(exists $seltx->{$ensid}){
		    $seltx->{$ensid}{$txid}=1;
		} else{
		    $seltx->{$ensid}={$txid=>1};
		}
	    }
	}
	
	$line=<$input>;
    }
}

##########################################################################

sub writeTranscriptGTF{
    my $gene=$_[0];
    my $tx=$_[1];
    my $txexons=$_[2];
    my $exoncoords=$_[3];
    my $source=$_[4];
    my $output=$_[5];
    
    my $chr="NA";
    my $strand="NA";
    
    my %hashcoords;
    
    foreach my $exon (keys %{$txexons->{$tx}}){
	my $thischr=$exoncoords->{$exon}{"chr"};
	my $thisstrand=$exoncoords->{$exon}{"strand"};
	my $thisstart=$exoncoords->{$exon}{"start"};
	my $thisend=$exoncoords->{$exon}{"end"};
	
	if($chr eq "NA"){
	    $chr=$thischr;
	} else{
	    if($chr ne $thischr){
		print "Weird! multiple chromosomes for ".$tx."\n";
		exit(1);
	    }
	}
	
	if($strand eq "NA"){
	    $strand=$thisstrand;
	} else{
	    if($strand ne $thisstrand){
		print "Weird! multiple strands for ".$tx."\n";
		exit(1);
	    }
	}

	if($thisend<$thisstart){
	    print "Weird coordinates for ".$exon."\n";
	    exit(1);
	}
	
	if(exists $hashcoords{$thisstart}){
	    print "Weird! already saw ".$thisstart." for ".$tx."\n";
	    exit(1);
	}
	
	$hashcoords{$thisstart}=$thisend;
    }
    
    my @sortedstart=sort {$a<=>$b} keys %hashcoords;

    foreach my $start (@sortedstart){
	my $end=$hashcoords{$start};

	print $output $chr."\t".$source."\texon\t".$start."\t".$end."\t.\t".$strand."\t.\tgene_id \"".$gene."\"; transcript_id \"".$tx."\";\n";
    }
}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script combines de novo and Ensembl annotations.\n";
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
$parameters{"pathEnsemblGTF"}="NA";
$parameters{"pathAssembledGTF"}="NA";
$parameters{"pathSelectedTranscripts"}="NA";
$parameters{"pathOutputGTF"}="NA";
$parameters{"pathOutputLog"}="NA";

my @defaultpars=("pathEnsemblGTF", "pathAssembledGTF", "pathSelectedTranscripts",  "pathOutputGTF", "pathOutputLog");


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

print "Reading Ensembl annotations...\n";

my %exoncoordsens;
my %exontxens;
my %genetxens;
my %txgeneens;
my %txexonens;

readGTF($parameters{"pathEnsemblGTF"}, \%exoncoordsens, \%exontxens, \%genetxens, \%txgeneens, \%txexonens);

my $nbexons=keys %exoncoordsens;
my $nbtx=keys %txgeneens;
my $nbg=keys %genetxens;

print "Found ".$nbexons." exons, ".$nbtx." transcripts and ".$nbg." genes in Ensembl annotations.\n";

##############################################################

print "Reading de novo annotations...\n";

my %exoncoordsass;
my %exontxass;
my %genetxass;
my %txgeneass;
my %txexonass;

readGTF($parameters{"pathAssembledGTF"}, \%exoncoordsass, \%exontxass, \%genetxass, \%txgeneass, \%txexonass);

my $nbexons=keys %exoncoordsass;
my $nbtx=keys %txgeneass;
my $nbg=keys %genetxass;

print "Found ".$nbexons." exons, ".$nbtx." transcripts and ".$nbg." genes in de novo annotations.\n";

##############################################################

print "Reading selected transcripts...\n";

my %selectedtx;
readSelectedTranscripts($parameters{"pathSelectedTranscripts"}, \%selectedtx);

my $nbgs=keys %selectedtx;
print "Found ".$nbgs." genes with selected transcripts.\n"; 

print "Done.\n";

##############################################################

print "Combining annotations and writing output...\n";

open(my $output, ">".$parameters{"pathOutputGTF"});
open(my $outlog, ">".$parameters{"pathOutputLog"});

foreach my $gene (keys %genetxens){
    foreach my $enstx (keys %{$genetxens{$gene}}){
	if(!(exists $selectedtx{$gene}{$enstx})){
	    writeTranscriptGTF($gene, $enstx, \%txexonens, \%exoncoordsens, "Ensembl", $output);
	} else{
	    print $outlog $enstx." ".$gene." was modified in StringTie annotations.\n";
	}
    }

    if(exists $selectedtx{$gene}){
	foreach my $asstx (keys %{$selectedtx{$gene}}){
	    if(exists $genetxens{$gene}{$asstx}){
		print $outlog "Also saw ".$asstx." in Ensembl annotations, keeping StringTie annotations.\n";
	    }
	    
	    writeTranscriptGTF($gene, $asstx, \%txexonass, \%exoncoordsass, "StringTie", $output);
	}
    }
}

## adding new genes

foreach my $gene (keys %selectedtx){
    if(!exists $genetxens{$gene}){
	foreach my $asstx (keys %{$selectedtx{$gene}}){
	    writeTranscriptGTF($gene, $asstx, \%txexonass, \%exoncoordsass, "StringTie", $output);
	}
    }
}

close($output);
close($outlog);

print "Done.\n";

##############################################################
