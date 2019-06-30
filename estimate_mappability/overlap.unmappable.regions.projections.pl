use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readProjectedExons{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $genesexons=$_[2];
    
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
	my $geneid=$s[$header{"GeneID"}];
		
	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0; ## 1-based
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	
	my $exonid=$chr.",".$start.",".$end.",".$strand;

	$exoncoords->{$exonid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	
	if(exists $genesexons->{$geneid}){
	    $genesexons->{$geneid}{$exonid}=1;
	}
	else{
	    $genesexons->{$geneid}={$exonid=>1};
	}
        
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub makeBlocks{
    my $genesexons=$_[0];
    my $exoncoords=$_[1];
    my $genesblocks=$_[2];
    my $blockcoords=$_[3];

    foreach my $gene (keys %{$genesexons}){
	my %hashpos;
	
	my $chr="NA";
	my $strand="NA";

	$genesblocks->{$gene}={};

	foreach my $exonid (keys %{$genesexons->{$gene}}){
	    my $start=$exoncoords->{$exonid}{"start"};
	    my $end=$exoncoords->{$exonid}{"end"};
	    
	    if($chr eq "NA"){
		$chr=$exoncoords->{$exonid}{"chr"};
		$strand=$exoncoords->{$exonid}{"strand"};
	    } else{
		if($exoncoords->{$exonid}{"chr"} ne $chr){
		    print "Saw two different chromosomes for ".$gene.": ".$chr." ".$exoncoords->{$exonid}{"chr"}."\n";
		    exit(1);
		}

		if($exoncoords->{$exonid}{"strand"} ne $strand){
		    print "Saw two different strands for ".$gene.": ".$strand." ".$exoncoords->{$exonid}{"strand"}."\n";
		    exit(1);
		}
	    }
	    
	    if(exists $hashpos{$start}){
		if($end>$hashpos{$start}){
		    $hashpos{$start}=$end;
		}
	    }
	    else{
		$hashpos{$start}=$end;
	    }
	}
    

	my @uniquepos=keys %hashpos;
	my @sortedpos=sort {$a<=>$b} @uniquepos;
	
	my $currentstart=$sortedpos[0];
	my $currentend=$hashpos{$sortedpos[0]};
	
	for(my $i=1; $i<@sortedpos; $i++){
	    my $thisstart=$sortedpos[$i];
	    my $thisend=$hashpos{$sortedpos[$i]};
	    
	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		if($currentend<$thisend){
		    $currentend=$thisend;
		}
	    }
	    else{
		my $currentid=$chr.",".$currentstart.",".$currentend.",".$strand;
		$blockcoords->{$currentid}={"chr"=>$chr, "start"=>$currentstart, "end"=>$currentend, "strand"=>$strand};
		$genesblocks->{$gene}{$currentid}=1;
		
		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}
	
	my $currentid=$chr.",".$currentstart.",".$currentend.",".$strand;
	$blockcoords->{$currentid}={"chr"=>$chr, "start"=>$currentstart, "end"=>$currentend, "strand"=>$strand};
	$genesblocks->{$gene}{$currentid}=1;
	
    }
}

##############################################################

sub readUnmappableRegions{
    my $pathin=$_[0];
    my $regions=$_[1];

    open(my $input, $pathin);
    ## no header
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[0];
	my $start=$s[1];
	my $end=$s[2];
	
	if(exists $regions->{$chr}){
	    push(@{$regions->{$chr}{"start"}}, $start);
	    push(@{$regions->{$chr}{"end"}}, $end);
	} else{
	    $regions->{$chr}={"start"=>[$start], "end"=>[$end]};
	}

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub orderExons{
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
	    
	    for(my $i=0; $i<$nbblocks; $i++){
		my $strand=${$hashstart{$chr}{$b}{"strand"}}[$i];
		
		push(@{$refordered->{$chr}{"start"}},$b);
		push(@{$refordered->{$chr}{"strand"}},$strand);
		push(@{$refordered->{$chr}{"end"}},${$hashstart{$chr}{$b}{"end"}}[$i]);
		push(@{$refordered->{$chr}{"id"}},${$hashstart{$chr}{$b}{"id"}}[$i]);
	    }
	}
    }
}

##############################################################

sub computeOverlap{
    my $orderedexons=$_[0];
    my $regions=$_[1];
    my $overlap=$_[2];
    
    foreach my $chr (keys %{$orderedexons}){
	if(exists $regions->{$chr}){
	    my $nbex=@{$orderedexons->{$chr}{"start"}};
	    my $nbreg=@{$regions->{$chr}{"start"}};
	    
	    my $firstreg=0;
	    
	    for(my $i=0; $i<$nbex; $i++){
		my $startex=${$orderedexons->{$chr}{"start"}}[$i];
		my $endex=${$orderedexons->{$chr}{"end"}}[$i];
		my $idex=${$orderedexons->{$chr}{"id"}}[$i];
	
		my $j=$firstreg;
		
		while($j<$nbreg && ${$regions->{$chr}{"end"}}[$j] < $startex){
		    $j++;
		}
		
		$firstreg=$j;
		
		while($j<$nbreg && ${$regions->{$chr}{"start"}}[$j] <= $endex){
		    my $startreg=${$regions->{$chr}{"start"}}[$j];
		    my $endreg=${$regions->{$chr}{"end"}}[$j];
		 		    
		    my $M=max($startex, $startreg);
		    my $m=min($endex, $endreg);
		    
		    if($M<=$m){
			if(exists $overlap->{$idex}){
			    push(@{$overlap->{$idex}{"start"}}, $M);
			    push(@{$overlap->{$idex}{"end"}}, $m);
			}
			else{
			    $overlap->{$idex}={"start"=>[$M], "end"=>[$m]};
			}
		    }
		    
		    $j++;
		}
	    }
	}
    }	
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes overlap with unmappable regions.\n";
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
$parameters{"pathProjectedExons"}="NA";
$parameters{"pathUnmappableRegions"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathProjectedExons", "pathUnmappableRegions", "pathOutput");

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

print "Reading projected annotations...\n";

my %exoncoords;
my %genesexons;

readProjectedExons($parameters{"pathProjectedExons"}, \%exoncoords, \%genesexons);

my $nbex=keys %exoncoords;
my $nbg=keys %genesexons;

print "Found ".$nbex." exons and ".$nbg." genes.\n";

print "Done.\n";

##############################################################

print "Making exon blocks...\n";

my %blockcoords;
my %genesblocks;

makeBlocks(\%genesexons, \%exoncoords, \%genesblocks, \%blockcoords);
  
print "Done.\n";

##############################################################

print "Reading unmappable regions..\n";

my %unmap;

readUnmappableRegions($parameters{"pathUnmappableRegions"}, \%unmap);

my $nbreg=keys %unmap;

print "Found unmappable regions on ".$nbreg." chromosomes.\n";

print "Done.\n";

##############################################################

print "Ordering coordinates...\n";

my %orderedblocks;

orderExons(\%blockcoords, \%orderedblocks);

print "Done.\n";

##############################################################

print "Computing overlap with unmappable regions...\n";

my %overlapunmap;

computeOverlap(\%orderedblocks, \%unmap,  \%overlapunmap);

my $nbov=keys %overlapunmap;

print "Found ".$nbov." exons that overlap with unmappable regions.\n";

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tTotalExonicLength\tUnmappableLength\n";

foreach my $geneid (keys %genesblocks){
    my $totlen=0;
    my $unmaplen=0;
    
    foreach my $blockid (keys %{$genesblocks{$geneid}}){
	my $start=$blockcoords{$blockid}{"start"};
	my $end=$blockcoords{$blockid}{"end"};
	
	$totlen+=($end-$start+1);
	
	if(exists $overlapunmap{$blockid}){
	   
	    my $nbex=@{$overlapunmap{$blockid}{"start"}};
	    
	    for(my $i=0; $i<$nbex; $i++){
		my $s=${$overlapunmap{$blockid}{"start"}}[$i];
		my $e=${$overlapunmap{$blockid}{"end"}}[$i];

		$unmaplen+=($e-$s+1);
	    }
	}
    }

    print $output $geneid."\t".$totlen."\t".$unmaplen."\n";
}

close($output);

print "Done.\n";

##############################################################
##############################################################
