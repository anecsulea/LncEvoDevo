use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

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

###########################################################################

sub extractTSS{
    my $transcriptcoords=$_[0];
    my $tsscoords=$_[1];
    
    foreach my $tx (keys %{$transcriptcoords}){
	my $chr=$transcriptcoords->{$tx}{"chr"};
	my $strand=$transcriptcoords->{$tx}{"strand"};
	my $start=$transcriptcoords->{$tx}{"start"};
	my $end=$transcriptcoords->{$tx}{"end"};
	my $gene=$transcriptcoords->{$tx}{"gene"};

	if($strand eq "1"){
	    $tsscoords->{$tx}={"gene"=>$gene, "chr"=>$chr, "strand"=>$strand, "position"=>$start};
	} else{
	    if($strand eq "-1"){
		$tsscoords->{$tx}={"gene"=>$gene, "chr"=>$chr, "strand"=>$strand, "position"=>$end};
	    }
	    else{
		print "Weird strand found for tx ".$tx.": ".$strand."\n";
		exit(1);
	    }
	}
    }
}

###########################################################################

sub orderTSS{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $tx (keys %{$coords}){
	my $chr=$coords->{$tx}{"chr"};
	my $strand=$coords->{$tx}{"strand"};
	my $position=$coords->{$tx}{"position"};
	
	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$position}){
		push(@{$hashpos{$chr}{$position}{"id"}},$tx);
		push(@{$hashpos{$chr}{$position}{"strand"}},$strand);
	    }
	    else{
		$hashpos{$chr}{$position}={"id"=>[$tx], "strand"=>[$strand]};
	    }
	}
	else{
	    $hashpos{$chr}={$position=>{"id"=>[$tx], "strand"=>[$strand]}};
	}
    }

    foreach my $chr (keys %hashpos){
	my @uniqueposition=keys %{$hashpos{$chr}};
	my @sortedposition=sort {$a<=>$b} @uniqueposition;

	foreach my $position (@sortedposition){
	    my $nb=@{$hashpos{$chr}{$position}{"id"}};

	    for(my $i=0; $i<$nb; $i++){
		my $strand=${$hashpos{$chr}{$position}{"strand"}}[$i];

		if(!(exists $ordered->{$chr})){
		    $ordered->{$chr}={"position"=>[], "strand"=>[], "id"=>[]};
		}

		my $id=${$hashpos{$chr}{$position}{"id"}}[$i];
	
		push(@{$ordered->{$chr}{"position"}}, $position);
		push(@{$ordered->{$chr}{"strand"}}, $strand);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

##############################################################################

sub computeDistanceAntisenseTSS{
    my $orderedcoords=$_[0]; 
    my $maxdist=$_[1];
    my $distance=$_[2];

    ## we compute the distance to the next exon on antisense

    foreach my $chr (keys %{$orderedcoords}){
	my $nb=@{$orderedcoords->{$chr}{"position"}};

	for(my $i=0; $i<$nb; $i++){
	    my $thisid=${$orderedcoords->{$chr}{"id"}}[$i];
	    my $thispos=${$orderedcoords->{$chr}{"position"}}[$i];
	    my $thisstrand=${$orderedcoords->{$chr}{"strand"}}[$i];
 
	    my $distleft="Inf";
	    my $idleft="NA";

	    my $j=$i-1;
	    
	    while($j>=0){
		my $nextid=${$orderedcoords->{$chr}{"id"}}[$j];
		my $nextstrand=${$orderedcoords->{$chr}{"strand"}}[$j];
		my $nextpos=${$orderedcoords->{$chr}{"position"}}[$j];
		
		$distleft=$thispos-$nextpos;
		
		if($nextstrand ne $thisstrand && $distleft <= $maxdist){
		    $idleft=$nextid;
		    last; ## we stop at the first opposite strand we see
		}

		if($distleft>$maxdist){
		    last;
		}
		
		$j--;
	    }
	    
	    my $distright="Inf";
	    my $idright="NA";
	    
	    my $j=$i+1;
	    
	    while($j<$nb){
		my $nextid=${$orderedcoords->{$chr}{"id"}}[$j];
		my $nextstrand=${$orderedcoords->{$chr}{"strand"}}[$j];
		my $nextpos=${$orderedcoords->{$chr}{"position"}}[$j];
		
		$distright=$nextpos-$thispos;
		
		if($nextstrand ne $thisstrand && $distright <= $maxdist){
		    $idright=$nextid;
		    last; ## we stop at the first opposite strand we see
		}
		
		if($distright>$maxdist){
		    last;
		}

		$j++;
	    }

	    if($idright eq "NA"){
		$distright=">".$maxdist;
	    }

	    if($idleft eq "NA"){
		$distleft=">".$maxdist;
	    }

	    $distance->{$thisid}={"distanceleft"=>$distleft, "closestleft"=>$idleft, "distanceright"=>$distright, "closestright"=>$idright};
	}
    }
}
    
###########################################################################

sub extractCloseAntisenseTSS{
    my $orderedcoords=$_[0]; 
    my $maxdistance=$_[1];
    my $closetss=$_[2];

    ## we compute the distance to the next exon on antisense

    foreach my $chr (keys %{$orderedcoords}){
	my $nb=@{$orderedcoords->{$chr}{"position"}};

	for(my $i=0; $i<$nb; $i++){
	    my $thisid=${$orderedcoords->{$chr}{"id"}}[$i];
	    my $thispos=${$orderedcoords->{$chr}{"position"}}[$i];
	    my $thisstrand=${$orderedcoords->{$chr}{"strand"}}[$i];
	    
	    $closetss->{$thisid}=[];

	    my $j=$i-1;
	    my $thisdist=0;
	    
	    while($j>=0 && $thisdist<$maxdistance){
		my $nextid=${$orderedcoords->{$chr}{"id"}}[$j];
		my $nextstrand=${$orderedcoords->{$chr}{"strand"}}[$j];
		my $nextpos=${$orderedcoords->{$chr}{"position"}}[$j];
		
		$thisdist=$thispos-$nextpos;		
		
		if($nextstrand ne $thisstrand){
		    if($thisdist<=$maxdistance){
			push(@{$closetss->{$thisid}}, $nextid);
		    } 
		}
		
		$j--;
	    }
	    
	
	    $j=$i+1;
	    $thisdist=0;
	    
	    while($j<$nb && $thisdist<$maxdistance){
		my $nextid=${$orderedcoords->{$chr}{"id"}}[$j];
		my $nextstrand=${$orderedcoords->{$chr}{"strand"}}[$j];
		my $nextpos=${$orderedcoords->{$chr}{"position"}}[$j];
		
		$thisdist=$nextpos-$thispos;

		if($nextstrand ne $thisstrand){
		    if($thisdist<=$maxdistance){
			push(@{$closetss->{$thisid}}, $nextid);
		    } 
		}
		
		$j++;
	    }
	}
    }
}

###########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes the distance to the closest antisense TSS. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

############################################################################
############################################################################

my %parameters;

$parameters{"pathGTF"}="NA";
$parameters{"maxDistance"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathGTF", "maxDistance", "pathOutput");

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

#####################################################################
#####################################################################

print "Reading coordinates for the first annotation set...\n";

my %exoncoords;
my %transcriptcoords;
my %genetx;

readGTF($parameters{"pathGTF"}, \%exoncoords, \%transcriptcoords, \%genetx);
  
my $nbgene=keys %genetx;
my $nbtx=keys %transcriptcoords;
my $nbex=keys %exoncoords;

print "Found ".$nbgene." genes, ".$nbtx." transcripts and ".$nbex." exons in the annotation.\n";

print "Done.\n";

#####################################################################

print "Extracting and ordering TSS coordinates...\n";

my %tss;
extractTSS(\%transcriptcoords, \%tss);

my %orderedtss;
orderTSS(\%tss, \%orderedtss);

print "Done.\n";

#####################################################################

print "Computing distance to the next antisense TSS...\n";


my $maxdist=$parameters{"maxDistance"}+0;

print "Maximum distance: ".$maxdist."\n";

my %distancetss;
computeDistanceAntisenseTSS(\%orderedtss, $maxdist,  \%distancetss);

print "Done.\n";

print "Extracting all close antisense TSS...\n";

my %closetss;

extractCloseAntisenseTSS(\%orderedtss, $maxdist, \%closetss);

print "Done.\n";

#####################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tTranscriptID\tChr\tStart\tEnd\tStrand\tTSS\tMinDistanceLeft\tClosestLeft\tMaxDistanceRight\tClosestRight\tTranscriptsCloseTSS\tGenesCloseTSS\n";

foreach my $gene (keys %genetx){
    foreach my $tx (keys %{$genetx{$gene}}){
	my $chr=$transcriptcoords{$tx}{"chr"};
	my $start=$transcriptcoords{$tx}{"start"};
	my $end=$transcriptcoords{$tx}{"end"};
	my $strand=$transcriptcoords{$tx}{"strand"};
	my $tss=$tss{$tx}{"position"};

	my $distleft=$distancetss{$tx}{"distanceleft"};
	my $closestleft=$distancetss{$tx}{"closestleft"};

	my $distright=$distancetss{$tx}{"distanceright"};
	my $closestright=$distancetss{$tx}{"closestright"};

	my $nbclose=@{$closetss{$tx}};

	my $txclose="NA";
	my $geneclose="NA";

	
	if($nbclose>0){
	    $txclose=join(",", @{$closetss{$tx}});
	    my %gc;

	    foreach my $t (@{$closetss{$tx}}){
		my $g=$transcriptcoords{$t}{"gene"};
		$gc{$g}=1;
	    }
	    
	    $geneclose=join(",",keys %gc);
	}

	print $output $gene."\t".$tx."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$tss."\t".$distleft."\t".$closestleft."\t".$distright."\t".$closestright."\t".$txclose."\t".$geneclose."\n";
    }
}

close($output);

print "Done.\n";

#####################################################################
