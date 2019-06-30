use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $genes=$_[1];
    
    open(my $input, $pathin);
    
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $geneid=$s[0];
	my $chr=$s[2];
	my $start=$s[3]+0; ## 1-based
	my $end=$s[4]+0;
	my $strand=$s[5];
	
	my $exonid=$chr.",".$start.",".$end.",".$strand;
	
	if(exists $genes->{$geneid}){
	    push(@{$genes->{$geneid}}, $exonid);
	}
	else{
	    $genes->{$geneid}=[$exonid];
	}
		
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readProjectedExons{
    my $pathin=$_[0];
    my $projectedexons=$_[1];
    
    open(my $input, $pathin);

    my $line=<$input>; ## header
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
	
	my $id=$s[$header{"ExonID"}];
	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];

	$projectedexons->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readSyntenyPredictions{
    my $pathin=$_[0];
    my $coords=$_[1];
    
    open(my $input, $pathin);

    my $line=<$input>; ## header
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
	
	my $id=$s[$header{"GeneID"}];
	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];

	$coords->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub checkRearrangedExons{
    my $projectedexons=$_[0]; ## one gene only
    my $filteredexons=$_[1];
    my $rejectedexons=$_[2];
    
    my $nbexons=keys %{$projectedexons};
    my @keyexons=keys %{$projectedexons};

    if($nbexons>1){
	my %rearrangedpairs;

	for(my $i=0; $i<($nbexons-1); $i++){
	    my $id1=$keyexons[$i];
	    my @s=split(",", $id1);
	    my $origstart1=$s[1]+0;
	    my $origend1=$s[2]+0;
	    my $origstrand1=$s[3];

	    my $projstart1=$projectedexons->{$id1}{"start"};
	    my $projend1=$projectedexons->{$id1}{"end"};
	    my $projstrand1=$projectedexons->{$id1}{"strand"};
	    
	    for(my $j=($i+1); $j<$nbexons; $j++){
		my $id2=$keyexons[$j];
		my @s=split(",", $id2);
		my $origstart2=$s[1]+0;
		my $origend2=$s[2]+0;

		my $projstart2=$projectedexons->{$id2}{"start"};
		my $projend2=$projectedexons->{$id2}{"end"};
		my $projstrand2=$projectedexons->{$id2}{"strand"};
		
		my $origsign=($origstart2-$origstart1);
		my $projsign=($projstart2-$projstart1);

		if($origsign!=0){
		    my $joinedsign=$origsign*$projsign;

		    if(($origstrand1 eq $projstrand1 && $joinedsign<0) || ($origstrand1 ne $projstrand1 && $joinedsign>0)){
			if(exists $rearrangedpairs{$id1}){
			    $rearrangedpairs{$id1}++;
			}
			else{
			    $rearrangedpairs{$id1}=1;
			}
			
			if(exists $rearrangedpairs{$id2}){
			    $rearrangedpairs{$id2}++;
			}
			else{
			    $rearrangedpairs{$id2}=1;
			}
		    }
		}
	    }
	}

	my $nbmultiple=0;
	my $nbone=0;
	my $nbrearranged=keys %rearrangedpairs;

	my $idmultiple="NA";

	foreach my $id (keys %rearrangedpairs){
	    my $nbpairs=$rearrangedpairs{$id};

	    if($nbpairs>1){
		$idmultiple=$id;
		$nbmultiple++;
	    }
	    else{
		$nbone++;
	    }
	}

	if($nbmultiple==1 && $nbone==($nbrearranged-1)){
	    ## we only rejected the exon that is in conflict with multiple exons

	     $rejectedexons->{$idmultiple}={};
	     
	     foreach my $key (keys %{$projectedexons->{$idmultiple}}){
		 $rejectedexons->{$idmultiple}{$key}=$projectedexons->{$idmultiple}{$key};
	     }
	     
	     foreach my $exonid (keys %{$projectedexons}){
		if($exonid ne $idmultiple){
		    $filteredexons->{$exonid}={};
		    
		    foreach my $key (keys %{$projectedexons->{$exonid}}){
			$filteredexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
		    }
		}
	    }
	}
	else{
	    ## we reject all conflicted exons
	     foreach my $exonid (keys %rearrangedpairs){
		 $rejectedexons->{$exonid}={};
		 
		 foreach my $key (keys %{$projectedexons->{$exonid}}){
		     $rejectedexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
		 }
	     }
	}

	## keep all exons that are in correct order with everything else (local rearrangements will be excluded)
	
	foreach my $exonid (keys %{$projectedexons}){
	    if(!(exists $rearrangedpairs{$exonid})){
		$filteredexons->{$exonid}={};
		
		foreach my $key (keys %{$projectedexons->{$exonid}}){
		    $filteredexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
		}
	    }
	}
    }
    else{
	foreach my $exonid (keys %{$projectedexons}){
	    $filteredexons->{$exonid}={};

	    foreach my $key (keys %{$projectedexons->{$exonid}}){
		$filteredexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
	    }
	}
    }
}

##############################################################

sub checkProjectedChromosomes{
    my $projectedexons=$_[0]; 
    my $filteredexons=$_[1];
    my $rejectedexons=$_[2];
    
    my %chromonb;

    foreach my $id (keys %{$projectedexons}){
	my $chr=$projectedexons->{$id}{"chr"};
	
	if(exists $chromonb{$chr}){
	    $chromonb{$chr}++;
	}
	else{
	    $chromonb{$chr}=1;
	}
    }

    my $nbprojected=keys %{$projectedexons};

    my $acceptedchr="NA";
    
    foreach my $chr (keys %chromonb){
	my $propchr=($chromonb{$chr}+0.0)/($nbprojected+0.0);

	if($propchr>0.5){
	    $acceptedchr=$chr;
	}
    }

    foreach my $exonid (keys %{$projectedexons}){
	my $chr=$projectedexons->{$exonid}{"chr"};
	
	if($chr eq $acceptedchr){
	    $filteredexons->{$exonid}={};
	    
	    foreach my $key (keys %{$projectedexons->{$exonid}}){
		$filteredexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
	    }
	}
	else{
	    $rejectedexons->{$exonid}={};
	    
	     foreach my $key (keys %{$projectedexons->{$exonid}}){
		$rejectedexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
	    }
	}
    }
}

##############################################################

sub checkProjectedStrands{
    my $projectedexons=$_[0]; 
    my $filteredexons=$_[1];
    my $rejectedexons=$_[2];
    
    my %strandnb;

    foreach my $id (keys %{$projectedexons}){
	my $strand=$projectedexons->{$id}{"strand"};
	
	if(exists $strandnb{$strand}){
	    $strandnb{$strand}++;
	}
	else{
	    $strandnb{$strand}=1;
	}
    }

    my $nbprojected=keys %{$projectedexons};

    my $acceptedstrand="NA";
    
    foreach my $strand (keys %strandnb){
	my $propstrand=($strandnb{$strand}+0.0)/($nbprojected+0.0);

	if($propstrand>0.5){
	    $acceptedstrand=$strand;
	}
    }

    foreach my $exonid (keys %{$projectedexons}){
	my $strand=$projectedexons->{$exonid}{"strand"};
	
	if($strand eq $acceptedstrand){
	    $filteredexons->{$exonid}={};
	    
	    foreach my $key (keys %{$projectedexons->{$exonid}}){
		$filteredexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
	    }
	}
	else{
	    $rejectedexons->{$exonid}={};
	    
	     foreach my $key (keys %{$projectedexons->{$exonid}}){
		$rejectedexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
	    }
	}
    }
}

##############################################################

sub checkIntronSize{
    my $projectedexons=$_[0]; ## one gene only
    my $maxratio=$_[1];
    my $maxaddedsize=$_[2];
    my $filteredexons=$_[3];
    my $rejectedexons=$_[4];
    my $verbose=$_[5];
    
    my $nbexons=keys %{$projectedexons};
    my @keyexons=keys %{$projectedexons};

    if($nbexons>1){
	my %wrongintronsize;

	for(my $i=0; $i<($nbexons-1); $i++){
	    my $id1=$keyexons[$i];
	    my @s=split(",", $id1);
	    my $origstart1=$s[1]+0;
	    my $origend1=$s[2]+0;
	    my $origstrand1=$s[3];

	    my $projstart1=$projectedexons->{$id1}{"start"};
	    my $projend1=$projectedexons->{$id1}{"end"};
	    my $projstrand1=$projectedexons->{$id1}{"strand"};
	    
	    for(my $j=($i+1); $j<$nbexons; $j++){
		my $id2=$keyexons[$j];
		my @s=split(",", $id2);
		my $origstart2=$s[1]+0;
		my $origend2=$s[2]+0;
		
		my $origintronsize=0;
		
		if($origend1 < $origstart2){
		    $origintronsize=$origstart2-$origend1-1;
		}
		
		if($origend2 < $origstart1){
		    $origintronsize=$origstart1-$origend2-1;
		}
		
		if($origintronsize>0){
		    my $projstart2=$projectedexons->{$id2}{"start"};
		    my $projend2=$projectedexons->{$id2}{"end"};
		    my $projstrand2=$projectedexons->{$id2}{"strand"};
		    
		    my $projintronsize=0;
		    
		    if($projend1<$projstart2){
			$projintronsize=$projstart2-$projend1-1;
		    }
		    
		    if($projend2<$projstart1){
			$projintronsize=$projstart1-$projend2-1;
		    }
		    
		    my $thisratio=$projintronsize/$origintronsize;
		    my $thisadded=$projintronsize-$origintronsize;

		    if(($projintronsize==0) || ($thisratio > $maxratio && $thisadded > $maxaddedsize)){ 

			if($verbose==1){
			    print "Found bad intron size for ".$id1." and ".$id2.": original size ".$origintronsize." new size ".$projintronsize." ratio ".$thisratio." added ".$thisadded."\n";
			    print $projstart1." ".$projend1." ".$projstart2." ".$projend2."\n";
			}
			
			if(exists $wrongintronsize{$id1}){
			    $wrongintronsize{$id1}++;
			}
			else{
			    $wrongintronsize{$id1}=1;
			}
			
			if(exists $wrongintronsize{$id2}){
			    $wrongintronsize{$id2}++;
			}
			else{
			    $wrongintronsize{$id2}=1;
			}
		    }
		}
	    }
	}
	
	
	my $nbmultiple=0;
	my $nbone=0;
	my $nbrearranged=keys %wrongintronsize;

	my $idmultiple="NA";

	foreach my $id (keys %wrongintronsize){
	    my $nbpairs=$wrongintronsize{$id};

	    if($nbpairs>1){
		$idmultiple=$id;
		$nbmultiple++;
	    }
	    else{
		$nbone++;
	    }
	}

	if($nbmultiple==1 && $nbone==($nbrearranged-1)){
	    ## we only rejected the exon that is in conflict with multiple exons

	     $rejectedexons->{$idmultiple}={};
	     
	     foreach my $key (keys %{$projectedexons->{$idmultiple}}){
		 $rejectedexons->{$idmultiple}{$key}=$projectedexons->{$idmultiple}{$key};
	     }
	     
	     foreach my $exonid (keys %{$projectedexons}){
		if($exonid ne $idmultiple){
		    $filteredexons->{$exonid}={};
		    
		    foreach my $key (keys %{$projectedexons->{$exonid}}){
			$filteredexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
		    }
		}
	    }
	}
	else{
	    ## we reject all conflicted exons
	     foreach my $exonid (keys %wrongintronsize){
		 $rejectedexons->{$exonid}={};
		 
		 foreach my $key (keys %{$projectedexons->{$exonid}}){
		     $rejectedexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
		 }
	     }
	}

	## keep all exons that are in correct relationship with everything else 
	
	foreach my $exonid (keys %{$projectedexons}){
	    if(!(exists $wrongintronsize{$exonid})){
		$filteredexons->{$exonid}={};
		
		foreach my $key (keys %{$projectedexons->{$exonid}}){
		    $filteredexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
		}
	    }
	}
    }
    else{
	foreach my $exonid (keys %{$projectedexons}){
	    $filteredexons->{$exonid}={};

	    foreach my $key (keys %{$projectedexons->{$exonid}}){
		$filteredexons->{$exonid}{$key}=$projectedexons->{$exonid}{$key};
	    }
	}
    }
}

##############################################################

sub readChromosomeCorrespondence{
    my $pathin=$_[0]; 
    my $corresp=$_[1]; 

    open(my $input, $pathin);

    my $line=<$input>; ## header
  
    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $ens=$s[0];
	my $ucsc=$s[1];

	$corresp->{$ens}=$ucsc;
    
	$line=<$input>;
    }
    
    close($input);
}


##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script filters projected genes.\n";
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
$parameters{"pathProjectedExons"}="NA";
$parameters{"pathSyntenyPredictions"}="NA";
$parameters{"syntenyRange"}="NA";
$parameters{"maxIntronSizeRatio"}="NA";
$parameters{"maxAddedIntronSize"}="NA";
$parameters{"pathOutputFilteredExons"}="NA";
$parameters{"pathOutputRejectedExons"}="NA";
$parameters{"pathOutputLog"}="NA";

my @defaultpars=("pathExonBlocks", "pathProjectedExons", "pathSyntenyPredictions", "syntenyRange", "maxIntronSizeRatio", "maxAddedIntronSize", "pathOutputFilteredExons", "pathOutputRejectedExons", "pathOutputLog");

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

print "Reading exon blocks...\n";

my %genes;

readExonBlocks($parameters{"pathExonBlocks"}, \%genes);

print "Done.\n";

##############################################################

print "Reading projected exons...\n";

my %projectedexons;

readProjectedExons($parameters{"pathProjectedExons"}, \%projectedexons);

print "Done.\n";

##############################################################

print "Reading synteny predictions...\n";

my %syntenypreds;

readSyntenyPredictions($parameters{"pathSyntenyPredictions"}, \%syntenypreds);

my $syntenyRange=$parameters{"syntenyRange"}+0;

print "We keep projected exons if they are within ".$syntenyRange." of the synteny-predicted coordinates.\n";

print "Done.\n";

##############################################################

print "Filtering and assigning projections...\n";

my $maxintronratio=$parameters{"maxIntronSizeRatio"}+0;
my $maxintronadded=$parameters{"maxAddedIntronSize"}+0;

print "Maximum intron size ratio: ".$maxintronratio."\n";
print "Maximum added intron size: ".$maxintronadded."\n";

open(my $outputok, ">".$parameters{"pathOutputFilteredExons"});
open(my $outputrejected, ">".$parameters{"pathOutputRejectedExons"});
open(my $outputlog, ">".$parameters{"pathOutputLog"});

print $outputok "GeneID\tExonID\tChr\tStart\tEnd\tStrand\n";
print $outputrejected "GeneID\tExonID\tChr\tStart\tEnd\tStrand\tRejectionReason\n";

foreach my $gene (keys %genes){
    my $chrsyn="NA";
    my $startsyn="NA";
    my $endsyn="NA";
    my $strandsyn="NA";
    my $hassyn=0;

    if(exists $syntenypreds{$gene}){
	$hassyn=1;
	$chrsyn=$syntenypreds{$gene}{"chr"};
	$startsyn=$syntenypreds{$gene}{"start"};
	$endsyn=$syntenypreds{$gene}{"end"};
	$strandsyn=$syntenypreds{$gene}{"strand"};
    }

    my %theseprojections;

    foreach my $exonid (@{$genes{$gene}}){
	if(exists $projectedexons{$exonid}){
	   
	    my $thischr=$projectedexons{$exonid}{"chr"};
	    my $thisstart=$projectedexons{$exonid}{"start"};
	    my $thisend=$projectedexons{$exonid}{"end"};
	    my $thisstrand=$projectedexons{$exonid}{"strand"};

	    my $M=max($startsyn-$syntenyRange, $thisstart);
	    my $m=min($endsyn+$syntenyRange, $thisend);

	    if(($hassyn==0) || ($hassyn==1 && $thischr eq $chrsyn && $thisstrand eq $strandsyn && $M<=$m)){
		$theseprojections{$exonid}={};
 
		foreach my $key (keys %{$projectedexons{$exonid}}){
		    $theseprojections{$exonid}{$key}=$projectedexons{$exonid}{$key};
		}
	    }
	    else{
		print $outputrejected $gene."\t".$exonid."\t".$thischr."\t".$thisstart."\t".$thisend."\t".$thisstrand."\tNotInSynteny\n";
	    }
	}
    }

    ## filter by chromosome

    my %filtered1;
    my %rejected1;

    checkProjectedChromosomes(\%theseprojections, \%filtered1, \%rejected1);
    
    ## filter by strand

    my %filtered2;
    my %rejected2;

    checkProjectedStrands(\%filtered1, \%filtered2, \%rejected2);

    ## filter rearrangements

    my %filtered3;
    my %rejected3;

    checkRearrangedExons(\%filtered2, \%filtered3, \%rejected3);


    ## filter intron size

    my %filtered4;
    my %rejected4;

    my $verbose=0;

    # test purposes
    # if($gene eq "ENSMUSG00000026774"){
    # 	$verbose=1;
    # }

    checkIntronSize(\%filtered3, $maxintronratio, $maxintronadded, \%filtered4, \%rejected4, $verbose);


    ## write output 

    foreach my $exonid (keys %filtered4){
	print $outputok $gene."\t".$exonid."\t".$filtered3{$exonid}{"chr"}."\t".$filtered3{$exonid}{"start"}."\t".$filtered3{$exonid}{"end"}."\t".$filtered3{$exonid}{"strand"}."\n";
    }

    foreach my $exonid (keys %rejected1){
	print $outputrejected $gene."\t".$exonid."\t".$rejected1{$exonid}{"chr"}."\t".$rejected1{$exonid}{"start"}."\t".$rejected1{$exonid}{"end"}."\t".$rejected1{$exonid}{"strand"}."\tWrongChromosome\n";
    }

    foreach my $exonid (keys %rejected2){
	print $outputrejected $gene."\t".$exonid."\t".$rejected2{$exonid}{"chr"}."\t".$rejected2{$exonid}{"start"}."\t".$rejected2{$exonid}{"end"}."\t".$rejected2{$exonid}{"strand"}."\tWrongStrand\n";
    }

    foreach my $exonid (keys %rejected3){
	print $outputrejected $gene."\t".$exonid."\t".$rejected3{$exonid}{"chr"}."\t".$rejected3{$exonid}{"start"}."\t".$rejected3{$exonid}{"end"}."\t".$rejected3{$exonid}{"strand"}."\tRearrangement\n";
    }

foreach my $exonid (keys %rejected4){
	print $outputrejected $gene."\t".$exonid."\t".$rejected4{$exonid}{"chr"}."\t".$rejected4{$exonid}{"start"}."\t".$rejected4{$exonid}{"end"}."\t".$rejected4{$exonid}{"strand"}."\tBadIntronSize\n";
    }
    
}

close($outputlog);
close($outputrejected);
close($outputok);

print "Done.\n";


##############################################################

## get consensus chromosome, start, end positions for each gene in each species
## flag wrongly assembled parts of chromosomes (genes split into multiple chromosomes)
## each projected exon gets a flag: gene; gene_wrong_assembly
## compare gene length original & projected: not too much length difference
## assign orthologous exons based on that



