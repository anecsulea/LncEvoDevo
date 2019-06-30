#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $transcriptcoords=$_[2];
    my $genetx=$_[3];
      
    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon"){
	    my $chr=$s[0];
	    my $start=$s[3]+0;
	    my $end=$s[4]+0;
	    my $ss=$s[6];
	    my $strand="NA";
	    
	    if($ss eq "+"){
		$strand="1";
	    } else{
		if($ss eq "-"){
		    $strand="-1";
		} 
	    }

	    ## we discard here unstranded transcripts
	    
	    if($strand ne "NA"){
		my $info=$s[8];
		my @infoarray=split(";", $info);
		
		my $gene=findInfo("gene_id", \@infoarray);
		my $tx=findInfo("transcript_id", \@infoarray);
		
		if($gene eq "NA"){
		    print "Weird! cannot find gene id for ".$line."\n";
		    exit(1);
		}
		
		if($tx eq "NA"){
		    print "Weird! cannot find transcript id for ".$line."\n";
		    exit(1);
		}
		
		if(exists $genetx->{$gene}){
		    $genetx->{$gene}{$tx}=1;
		} else{
		    $genetx->{$gene}={$tx=>1};
		}

		my $exonid=$chr.",".$start.",".$end.",".$strand;
		
		$exoncoords->{$exonid}={"start"=>$start, "end"=>$end, "chr"=>$chr, "strand"=>$strand};
		
		if(exists $transcriptcoords->{$tx}){
		    if($start < $transcriptcoords->{$tx}{"start"}){
			$transcriptcoords->{$tx}{"start"}=$start;
		    }
		    
		    if($end > $transcriptcoords->{$tx}{"end"}){
			$transcriptcoords->{$tx}{"end"}=$end;
		    }
		    
		} else{
		    $transcriptcoords->{$tx}={"start"=>$start, "end"=>$end, "chr"=>$chr, "strand"=>$strand, "gene"=>$gene};
		}
	    }
	}
	
	$line=<$input>;
    }
    
    close($input);
}

#########################################################################

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

##########################################################################

sub extractPromoters{
    my $genetx=$_[0];
    my $transcriptcoords=$_[1];
    my $size=$_[2];
    my $promoters=$_[3];

    ## we take only one promoter per gene, the most distal one

    foreach my $gene (keys %{$genetx}){
	foreach my $tx (keys %{$genetx->{$gene}}){
	    my $chr=$transcriptcoords->{$tx}{"chr"};
	    my $strand=$transcriptcoords->{$tx}{"strand"};
	    
	    my $start=$transcriptcoords->{$tx}{"start"};
	    my $end=$transcriptcoords->{$tx}{"end"};
	    
	    if($strand eq "1" || $strand eq "+"){
		if(exists $promoters->{$gene}){
		    my $prevstart=${$promoters->{$gene}{"end"}}[0];
		    
		    if($start<$prevstart){
			$promoters->{$gene}={"chr"=>$chr, "start"=>[$start-$size], "end"=>[$start], "strand"=>$strand};
		    }
		}
		else{
		    $promoters->{$gene}={"chr"=>$chr, "start"=>[$start-$size], "end"=>[$start], "strand"=>$strand};
		}
		
	    }
	    
	    if($strand eq "-1" || $strand eq "+"){
		if(exists $promoters->{$gene}){
		    my $prevstart=${$promoters->{$gene}{"start"}}[0];

		    if($end>$prevstart){ 
			$promoters->{$gene}={"chr"=>$chr, "start"=>[$end], "end"=>[$end+$size], "strand"=>$strand};
		    }
		} else{
		    $promoters->{$gene}={"chr"=>$chr, "start"=>[$end], "end"=>[$end+$size], "strand"=>$strand};
		}
	    }
	}
    }
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
##########################################################################

sub overlapBlocks{
    my $refblocks1=$_[0]; ## ordered exon blocks
    my $refblocks2=$_[1]; ## ordered exon blocks
    my $refoverlap=$_[2]; ## overlap with blocks
   
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
		    
		    my $classfam=${$refblocks2->{$chr}{"classfam"}}[$j];
		    my $name=${$refblocks2->{$chr}{"name"}}[$j];
		    my $strand2=${$refblocks2->{$chr}{"strand"}}[$j];

		    if($classfam eq ""){
			print "Found null type repeat : ".$chr." ".${$refblocks2->{$chr}{"start"}}[$j]." ".${$refblocks2->{$chr}{"end"}}[$j]."\n";
			exit(1);
		    }
		    
		    my $lenov=$m-$M+1;
		    
		    if($lenov>=1){
			
			my $keyblock=$start1.",".$end1;
			
			my $keyblock2=${$refblocks2->{$chr}{"start"}}[$j]." ".${$refblocks2->{$chr}{"end"}}[$j];
			
			my $frov1=($m-$M+1.0)/($end1-$start1+1.0);
			my $frov2=($m-$M+1.0)/(${$refblocks2->{$chr}{"end"}}[$j]-${$refblocks2->{$chr}{"start"}}[$j]+1.0);
			
			if(exists $refoverlap->{$gene1}){
			    if(!(exists $refoverlap->{$gene1}{$keyblock})){
				$refoverlap->{$gene1}{$keyblock}={"start"=>[$M],"end"=>[$m],"name"=>[$name],"classfam"=>[$classfam]};
			    }
			    else{
				push(@{$refoverlap->{$gene1}{$keyblock}{"start"}},$M);
				push(@{$refoverlap->{$gene1}{$keyblock}{"end"}},$m);
				push(@{$refoverlap->{$gene1}{$keyblock}{"classfam"}},$classfam);
				push(@{$refoverlap->{$gene1}{$keyblock}{"name"}},$name);
			    }
			}
			else{
			    $refoverlap->{$gene1}={$keyblock=>{"start"=>[$M],"end"=>[$m],"name"=>[$name],"classfam"=>[$classfam]}};
			}
		    }
		    
		    
		    $j++;
		}
	    }
	}
    }
}

##########################################################################

sub readRepeatMasker{
    my $pathin=$_[0];
    my $type=$_[1];
    my $okclasses=$_[2];
    my $repeatmasker=$_[3];

    my %hashrepeats;
    my @s=split("\\.", $pathin);
    my $ext=$s[-1];
    my $input;

    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    } else{
	open($input,$pathin);
    }

    
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    my $indexclass="NA";

    if(exists $header{$type}){
	$indexclass=$header{$type};
    } else{
	print "Cannot find column corresponding to ".$type."\n";
	exit(1);
    }
 
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split(" ",$line);
	
	my $chr=$s[5];
	
	my $prefix=substr $chr,0,3;

	if($prefix eq "chr"){
	    $chr=substr $chr,3;
	}

	
	my $start=$s[6]+1; ## coordinates start at 1 now
	my $end=$s[7]+0;
	my $strand=$s[9];

	if($strand eq "+"){
	    $strand="1";
	}
	else{
	    $strand="-1";
	}
	
	my $name=$s[10];
	my $classfam=$s[$indexclass];

	if($classfam eq ""){
	    print $line."\n";
	    exit(1);
	}
	
	if(exists $okclasses->{$classfam} || exists $okclasses->{"any"}){
	    if(exists $hashrepeats{$chr}){
		if(exists $hashrepeats{$chr}{$start}){
		    push(@{$hashrepeats{$chr}{$start}{"end"}},$end);
		    push(@{$hashrepeats{$chr}{$start}{"name"}},$name);
		    push(@{$hashrepeats{$chr}{$start}{"classfam"}},$classfam);
		    push(@{$hashrepeats{$chr}{$start}{"strand"}},$strand);
		}
		else{
		    $hashrepeats{$chr}{$start}={"end"=>[$end],"name"=>[$name],"classfam"=>[$classfam],"strand"=>[$strand]};
		}
	    }
	    else{
		$hashrepeats{$chr}={$start=>{"end"=>[$end],"name"=>[$name],"classfam"=>[$classfam],"strand"=>[$strand]}};
	    }
	    
	}
	
	$line=<$input>;
    }
    
    close($input);

    ### now order repeats
    
    foreach my $chr (keys %hashrepeats){
	$repeatmasker->{$chr}={"start"=>[],"end"=>[],"name"=>[],"classfam"=>[],"strand"=>[]};
	
	my @uniquestart=keys %{$hashrepeats{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashrepeats{$chr}{$start}{"end"}};
	    for(my $i=0; $i<$nb; $i++){

		if(${$hashrepeats{$chr}{$start}{"classfam"}}[$i] eq ""){
		    print "Found null type repeat while sorting: ".$chr." ".$start." ".${$hashrepeats{$chr}{$start}{"end"}}[$i]."\n";
		    exit(1);
		}

		push(@{$repeatmasker->{$chr}{"start"}},$start);
		push(@{$repeatmasker->{$chr}{"end"}},${$hashrepeats{$chr}{$start}{"end"}}[$i]);
		push(@{$repeatmasker->{$chr}{"name"}},${$hashrepeats{$chr}{$start}{"name"}}[$i]);
		push(@{$repeatmasker->{$chr}{"classfam"}},${$hashrepeats{$chr}{$start}{"classfam"}}[$i]);
		push(@{$repeatmasker->{$chr}{"strand"}},${$hashrepeats{$chr}{$start}{"strand"}}[$i]);
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
    print "This script extracts overlaps between promoters and transposable elements \n";
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
$parameters{"promoterSize"}="NA";
$parameters{"pathRepeatMasker"}="NA";
$parameters{"type"}="NA";
$parameters{"repeatClasses"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathGTF","promoterSize", "pathRepeatMasker", "type", "repeatClasses","pathOutput");

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

print "Reading annotations...\n";

my %exoncoords;
my %transcriptcoords;
my %genetx;

readGTF($parameters{"pathGTF"}, \%exoncoords, \%transcriptcoords, \%genetx);

my $nbex=keys %exoncoords;
my $nbtx=keys %transcriptcoords;
my $nbg=keys %genetx;

print "Found ".$nbex." exons, ".$nbtx." transcripts and ".$nbg." genes.\n";
print "Done.\n";

#############################################################################

my $promsize=$parameters{"promoterSize"}+0;

print "Promoter size: ".$promsize."\n";

print "Extracting promoters...\n";
my %promoters;

extractPromoters(\%genetx, \%transcriptcoords, $promsize, \%promoters);

print "Done.\n";

print "Ordering promoters...\n";
my %orderedpromoters;
orderExonBlocks(\%promoters,\%orderedpromoters);
print "Done.\n";

##########################################################################

print "Reading RepeatMasker...\n";

my %classes;

my @cl=split(",",$parameters{"repeatClasses"});

foreach my $c (@cl){
    $classes{$c}=1;
}

print "Analyzing repeat classes: ".join(", ",keys %classes)."\n";

my %repeats;
my $type=$parameters{"type"};

readRepeatMasker($parameters{"pathRepeatMasker"}, $type, \%classes,\%repeats);

print "Done.\n";

##########################################################################

print "Extracting overlaps...\n";

my %overlap;
overlapBlocks(\%orderedpromoters,\%repeats,\%overlap);

print "Done.\n";

### we now sort the overlapping regions according to the type of the transposable element

my %classoverlap;

foreach my $gene (keys %overlap){
    foreach my $keyblock (keys %{$overlap{$gene}}){
	my $nbov=@{$overlap{$gene}{$keyblock}{"start"}};

	for(my $i=0;$i<$nbov;$i++){
	    my $start=${$overlap{$gene}{$keyblock}{"start"}}[$i];
	    my $end=${$overlap{$gene}{$keyblock}{"end"}}[$i];
	    my $name=${$overlap{$gene}{$keyblock}{"name"}}[$i];
	    my $class=${$overlap{$gene}{$keyblock}{"classfam"}}[$i];
	    
	    ## class

	     if(exists $classoverlap{$class}){
		if(exists $classoverlap{$class}{$gene}){
		    if(exists $classoverlap{$class}{$gene}{$keyblock}){
			push(@{$classoverlap{$class}{$gene}{$keyblock}{"start"}},$start);
			push(@{$classoverlap{$class}{$gene}{$keyblock}{"end"}},$end);
			push(@{$classoverlap{$class}{$gene}{$keyblock}{"name"}},$name);
		
		    }
		    else{
			$classoverlap{$class}{$gene}{$keyblock}={"start"=>[$start],"end"=>[$end],"name"=>[$name]};
		    }
		}
		else{
		    $classoverlap{$class}{$gene}={$keyblock=>{"start"=>[$start],"end"=>[$end],"name"=>[$name]}};
		}
	    }
	    else{
		$classoverlap{$class}={$gene=>{$keyblock=>{"start"=>[$start],"end"=>[$end],"name"=>[$name]}}};
	    }

	}
    }
}

##########################################################################

print "Making overlap blocks...\n";

my %classblocks;

foreach my $cl (keys %classoverlap){
    $classblocks{$cl}={};
    makeOverlapBlocks($classoverlap{$cl},  $classblocks{$cl});
}

print "Done.\n";

##########################################################################

print "Writing output...\n";

my @classes=keys %classblocks;
my @sortedclasses = sort @classes;


open(my $output,">".$parameters{"pathOutput"});

my $line="GeneID";
foreach my $fam (@sortedclasses){
    $line.="\t".$fam;
}

print $output $line."\n";

foreach my $gene (keys %promoters){
    my $chr=$promoters{$gene}{"chr"};
    my $nbpromoters=@{$promoters{$gene}{"start"}};
    
    my $totlength=0;

    for(my $i=0;$i<$nbpromoters;$i++){
	$totlength+=(${$promoters{$gene}{"end"}}[$i]-${$promoters{$gene}{"start"}}[$i]+1);
    }

    my $line=$gene;

    foreach my $fam (@sortedclasses){
	
	my $overlaplength=0;
    
	if(exists $classblocks{$fam}{$gene}){
	    
	    foreach my $idpromoter (keys %{$classblocks{$fam}{$gene}}){
		my $nbpieces=@{$classblocks{$fam}{$gene}{$idpromoter}{"start"}};
	    
		for(my $i=0;$i<$nbpieces;$i++){
		    $overlaplength+=(${$classblocks{$fam}{$gene}{$idpromoter}{"end"}}[$i]-${$classblocks{$fam}{$gene}{$idpromoter}{"start"}}[$i]+1);
		}
	    }
	}

	my $propoverlap=($overlaplength+0.0)/($totlength+0.0);

	$line.="\t".$propoverlap;
    }

    print $output $line."\n";

}
close($output);


print "Done.\n";

##########################################################################
