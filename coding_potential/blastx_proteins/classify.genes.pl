#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readGTF{
    my $pathin=$_[0];
    my $transcripts=$_[1];
    my $genestx=$_[2];
    
    open(my $input, $pathin);
    
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon"){
	    my $chr=$s[0];
	    my $start=$s[3]+0; ## 1-based
	    my $end=$s[4]+0;
	    my $strand=$s[6];
	  
	    my $newstrand=$strand;

	    if($strand eq "+"){
		$newstrand="1";
	    }
	    
	    if($strand eq "-"){
		$newstrand="-1";
	    }
	    	    
	    my $exonid=$chr.",".$start.",".$end.",".$newstrand;

	    my $info=$s[8];
	    my @t=split(";", $info);
	    my $txid=findInfo("transcript_id", \@t);
	    my $geneid=findInfo("gene_id", \@t);

	    if(exists $genestx->{$geneid}){
		$genestx->{$geneid}{$txid}=1;	
	    } else{
		$genestx->{$geneid}={$txid=>1};	
	    }

	    if(exists $transcripts->{$txid}){
		push(@{$transcripts->{$txid}{"start"}}, $start);
		push(@{$transcripts->{$txid}{"end"}}, $end);
	    }
	    else{
		$transcripts->{$txid}={"chr"=>$chr, "strand"=>$newstrand, "start"=>[$start], "end"=>[$end]};
	    }
	}

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];
    
    my $res="NA";
    my $nbfound=0;
    
    my @grepres=grep(/${pattern}/,@{$array});
    
    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
	$nbfound++;
    }
    else{
	foreach my $g (@grepres){
	    my @u=split(" ",$g);
	    
	    if($u[0] eq $pattern){
		$nbfound++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }
    
    if($nbfound==1){
	return $res;
    } else{
	return "NA";
    }
}

##############################################################

sub readBlastX{
    my $pathin=$_[0];
    my $maxpval=$_[1];
    my $minpcid=$_[2];
    my $blastx=$_[3];

    open(my $input, $pathin);
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $txid=$s[0];
	my $pcid=$s[2]+0.0;
	my $pvalue=$s[10]+0.0;

	my $startquery=$s[6]+0;
	my $endquery=$s[7]+0;

	if($pvalue<=$maxpval && $pcid>=$minpcid){
	    if(exists $blastx->{$txid}){
		push(@{$blastx->{$txid}{"start"}}, $startquery);
		push(@{$blastx->{$txid}{"end"}}, $endquery);
	    } else{
		$blastx->{$txid}={"start"=>[$startquery], "end"=>[$endquery]};
	    }
	}
	
	$line=<$input>;
    }

    close($input);
}


##############################################################

sub makeBlastBlocks{
    my $blasthits=$_[0];
    my $blocks=$_[1];

    foreach my $txid (keys %{$blasthits}){
	$blocks->{$txid}={"start"=>[], "end"=>[]};

	my $nbhits=@{$blasthits->{$txid}{"start"}};
	my %hashpos;
	
	for(my $i=0; $i<$nbhits; $i++){
	    my $start=${$blasthits->{$txid}{"start"}}[$i];
	    my $end=${$blasthits->{$txid}{"end"}}[$i];

	    if(exists $hashpos{$start}){
		if($end>$hashpos{$start}){
		    $hashpos{$start}=$end;
		}
	    } else{
		$hashpos{$start}=$end;
	    }
	}

	my @sortedstart=sort {$a<=>$b} (keys %hashpos);

	my $currentstart=$sortedstart[0];
	my $currentend=$hashpos{$currentstart};

	my $nbsorted=@sortedstart;

	if($nbsorted>0){
	    for(my $i=1; $i<@sortedstart; $i++){
		my $thisstart=$sortedstart[$i];
		my $thisend=$hashpos{$thisstart};
		
		if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		    if($thisend>$currentend){
			$currentend=$thisend;
		    }
		} else{
		    push(@{$blocks->{$txid}{"start"}}, $currentstart);
		    push(@{$blocks->{$txid}{"end"}}, $currentend);
		    
		    $currentstart=$thisstart;
		    $currentend=$thisend;
		}
	    }
	    
	    ## last block
	    
	    push(@{$blocks->{$txid}{"start"}}, $currentstart);
	    push(@{$blocks->{$txid}{"end"}}, $currentend);
	}
    }

}

##############################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
   
    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];

	    # print "saw chromosome ".$id."\n";
	    
	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script classifies genes as coding or noncoding based on blastx results.\n";
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

$parameters{"pathGTF"}="NA";
$parameters{"pathFastacDNAs"}="NA";
$parameters{"pathBlastX"}="NA";
$parameters{"maxEValue"}="NA";
$parameters{"minPCIdentity"}="NA";
$parameters{"minFractionOverlap"}="NA";
$parameters{"minLengthOverlap"}="NA";
$parameters{"pathOutputGenes"}="NA";
$parameters{"pathOutputTranscripts"}="NA";

my %defaultvalues;
my @defaultpars=("pathGTF", "pathFastacDNAs","pathBlastX", "maxEValue", "minPCIdentity", "minFractionOverlap", "minLengthOverlap", "pathOutputGenes", "pathOutputTranscripts");

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

#####################################################################################
#####################################################################################

print "Reading annotations...\n";

my %transcripts;
my %genes;

readGTF($parameters{"pathGTF"}, \%transcripts, \%genes);

my $nbtx=keys %transcripts;
my $nbg=keys %genes;

print "Found ".$nbtx." transcripts and ".$nbg." genes.\n";

print "Done.\n";

######################################################################################

print "Reading cDNA sequences...\n";

my %txseq;
readFasta($parameters{"pathFastacDNAs"}, \%txseq);
my $nbt=keys %txseq;

print "Found sequences for ".$nbt." transcripts.\n";

my %txlen;
foreach my $id (keys %txseq){
    my $seq=$txseq{$id};
    my $tot=length $seq;
    my $nbn = ($seq =~ tr/N//);
    my $unmasked=$tot-$nbn;

    $txlen{$id}={"total"=>$tot, "unmasked"=>$unmasked};
}
print "Done.\n";

######################################################################################

print "Reading blastx results...\n";

my $minpcid=$parameters{"minPCIdentity"}+0.0;
my $maxeval=$parameters{"maxEValue"}+0.0;

print "Minimum % identity: ".$minpcid.".\n";
print "Maximum e-value: ".$maxeval."\n";

my %blastx;

readBlastX($parameters{"pathBlastX"}, $maxeval, $minpcid, \%blastx); 

my $nbres=keys %blastx;

print "Found ".$nbres." transcripts with blastx hits.\n";

print "Making blast blocks...\n";

my %blastblocks;

makeBlastBlocks(\%blastx, \%blastblocks);

print "Done.\n";

#####################################################################################

print "Classifying genes and writing output...\n";

open(my $outputtx, ">".$parameters{"pathOutputTranscripts"});
open(my $outputgenes, ">".$parameters{"pathOutputGenes"});

my $minfroverlap=$parameters{"minFractionOverlap"}+0.0;
my $minlenoverlap=$parameters{"minLengthOverlap"}+0;

print "Minimum length overlap: ".$minlenoverlap." minimum fraction overlap: ".$minfroverlap."\n";

print $outputgenes "GeneID\tClass\tCodingTranscripts\n";
print $outputtx "GeneID\tTranscriptID\tTotalLength\tUnmaskedLength\tLengthBlastXOverlap\tClass\n";

foreach my $gene (keys %genes){
    my $classgene="NA";
    my %codingtranscripts;
    my %noncodingtranscripts;
    
    foreach my $tx (keys %{$genes{$gene}}){
	my $classtx="NA";
	
	if(!exists $txlen{$tx}){
	    print "Weird! cannot find transcript ".$tx." for ".$gene."\n";
	    exit(1);
	}
	
	my $totlen=$txlen{$tx}{"total"};
	my $unmaskedlen=$txlen{$tx}{"unmasked"};
	
	my $lenpositive=0;
	
	if(exists $blastblocks{$tx}){
	    my $nbb=@{$blastblocks{$tx}{"start"}};
	    for(my $i=0; $i<$nbb; $i++){
		$lenpositive+=(${$blastblocks{$tx}{"end"}}[$i]-${$blastblocks{$tx}{"start"}}[$i]+1);
	    }
	}
	
	if($totlen==0){
	    print "Weird! null length for ".$tx." and ".$gene."\n";
	    exit(1);
	}
	
	if($unmaskedlen>=$minlenoverlap){
	    my $frpos=($lenpositive+0.0)/($unmaskedlen+0.0);
	    
	    if($frpos>=$minfroverlap && $lenpositive>=$minlenoverlap){
		$classtx="coding";
		$classgene="coding";
		$codingtranscripts{$tx}=1;
	    } else{
		$classtx="noncoding";
		$noncodingtranscripts{$tx}=1;
	    }
	    
	    print $outputtx $gene."\t".$tx."\t".$totlen."\t".$unmaskedlen."\t".$lenpositive."\t".$classtx."\n";
	}
    }
    
    my $idcodingtx="NA";
    
    my $nbcoding=keys %codingtranscripts;
    my $nbnoncoding=keys %noncodingtranscripts;
    
    if($nbcoding>0){
	$idcodingtx=join(",",sort(keys %codingtranscripts));
    }
    else{
	if($nbnoncoding>0){
	    $classgene="noncoding";
	}
    }

    print $outputgenes $gene."\t".$classgene."\t".$idcodingtx."\n";
}

close($outputtx);
close($outputgenes);

print "Done.\n";

#####################################################################################
