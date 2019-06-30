use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $geneexons=$_[2];
 
    open(my $input, $pathin);
    
    my $line=<$input>;
    my $prefix=substr $line, 0,1;
    
    while($prefix eq "#"){
	$line=<$input>;
	$prefix=substr $line, 0,1;
    }

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $geneid=$s[0];
	my $chr=$s[2];
	my $start=$s[3]+0; ## 1-based
	my $end=$s[4]+0;
	my $strand=$s[5];
				
	my $exonid=$geneid.",".$chr.",".$start.",".$end.",".$strand;
	
	## fill in exon coords
	
	$exoncoords->{$exonid}={"chr"=>$chr,"start"=>$start, "end"=>$end, "strand"=>$strand, "gene"=>$geneid};
	
	if(exists $geneexons->{$geneid}){
	    push(@{$geneexons->{$geneid}}, $exonid);
	}
	else{
	    $geneexons->{$geneid}=[$exonid];
	}

	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub extractGeneCoords{
    my $geneexons=$_[0];
    my $exoncoords=$_[1];
    my $genecoords=$_[2];

    foreach my $gene (keys %{$geneexons}){
	my $chr="NA";
	my $start="NA";
	my $end="NA";
	my $strand="NA";
	
	foreach my $exon (@{$geneexons->{$gene}}){
	    my $thischr=$exoncoords->{$exon}{"chr"};
	    my $thisstart=$exoncoords->{$exon}{"start"};
	    my $thisend=$exoncoords->{$exon}{"end"};
	    my $thisstrand=$exoncoords->{$exon}{"strand"};

	    if($chr eq "NA"){
		$chr=$thischr;
	    } else{
		if($chr ne $thischr){
		    print "Weird! multiple chromosomes seen for ".$gene."\n";
		    exit(1);
		}
	    }

	    if($strand eq "NA"){
		$strand=$thisstrand;
	    } else{
		if($strand ne $thisstrand){
		    print "Weird! multiple strands seen for ".$gene.": ".$strand." ".$thisstrand."\n";
		    exit(1);
		}
	    }

	    if($start eq "NA"){
		$start=$thisstart;
	    } else{
		if($thisstart<$start){
		    $start=$thisstart;
		}
	    }

	    if($end eq "NA"){
		$end=$thisend;
	    } else{
		if($thisend>$end){
		    $end=$thisend;
		}
	    }
	}

	$genecoords->{$gene}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
    }
}

##############################################################

sub makeExonBlocks{
    my $genes=$_[0];
    my $exons=$_[1];
    my $collapse=$_[2];
    my $exonblocks=$_[3];

    foreach my $gene (keys %{$genes}){
	my $chr="NA";
	my $strand="NA";
	
	my %hashstart;
		
	foreach my $exon (@{$genes->{$gene}}){
	    if($chr eq "NA"){
		$chr=$exons->{$exon}{"chr"};
		$strand=$exons->{$exon}{"strand"};
	    }
	    else{
		if($exons->{$exon}{"chr"} ne $chr || $exons->{$exon}{"strand"} ne $strand){
		    print "Weird!!! saw chromosomes ".$exons->{$exon}{"chr"}." and ". $chr.", strands ". $exons->{$exon}{"strand"}." and ". $strand." for ".$gene."\n";
		    exit(1);
		}
	    }
	    
	    my $b=$exons->{$exon}{"start"};
	    my $e=$exons->{$exon}{"end"};
	    
	    if(exists $hashstart{$b}){
		$hashstart{$b}{$e}=1;
	    }
	    else{
		$hashstart{$b}={$e=>1};
	    }
	}	    
	  

	$exonblocks->{$gene}={"chr"=>$chr,"strand"=>$strand,"start"=>[],"end"=>[]};
	
	### now make exon blocks

	my @uniquestart = keys %hashstart;
	my @sortedstart = sort {$a <=> $b} @uniquestart;
	
	my @start;
	my @end;
	
	foreach my $beg (@sortedstart){
	    my @thisend=keys %{$hashstart{$beg}};
	    my @sortedend = sort {$a <=> $b} @thisend;
	    
	    foreach my $en (@sortedend){
		push(@start, $beg);
		push(@end, $en);
	    }
	}	
	
	my $nbex=@start;
	
	my $currentstart=$start[0];
	my $currentend=$end[0];
	
	for(my $u=1;$u<$nbex;$u++){
	    my $thisstart=$start[$u];
	    my $thisend=$end[$u];
	    		
	    ## cluster blocks if they overlap
	    
	    if($thisstart>=$currentstart && $thisstart<=($currentend+$collapse)){  
		
		## we only change the end if it's larger than the current position
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    }
	    else{
		push(@{$exonblocks->{$gene}{"start"}},$currentstart);
		push(@{$exonblocks->{$gene}{"end"}},$currentend);
		
		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}
	
	## don't forget the last block
	
	push(@{$exonblocks->{$gene}{"start"}},$currentstart);
	push(@{$exonblocks->{$gene}{"end"}},$currentend);
    }
}

##############################################################

sub orderExons{
    my $exons=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $exid (keys %{$exons}){
	my $chr=$exons->{$exid}{"chr"};
	my $strand=$exons->{$exid}{"strand"};
	my $start=$exons->{$exid}{"start"};
	my $end=$exons->{$exid}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		push(@{$hashpos{$chr}{$start}{"id"}},$exid);
		push(@{$hashpos{$chr}{$start}{"end"}},$end);
		push(@{$hashpos{$chr}{$start}{"strand"}},$strand);
	    }
	    else{
		$hashpos{$chr}{$start}={"id"=>[$exid],"end"=>[$end],"strand"=>[$strand]};
	    }
	}
	else{
	    $hashpos{$chr}={$start=>{"id"=>[$exid],"end"=>[$end],"strand"=>[$strand]}};
	}
    }

    foreach my $chr (keys %hashpos){
	$ordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"id"=>[]};

	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashpos{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		my $end=${$hashpos{$chr}{$start}{"end"}}[$i];
		my $strand=${$hashpos{$chr}{$start}{"strand"}}[$i];
		my $id=${$hashpos{$chr}{$start}{"id"}}[$i];

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"strand"}}, $strand);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

##############################################################

sub extractOverlap{
    my $coords1=$_[0]; ## ordered coordinates
    my $coords2=$_[1]; ## ordered coordinates
    my $margin=$_[2];
    my $type=$_[3];
    my $overlap=$_[4];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nbex1=@{$coords1->{$chr}{"start"}};
	    my $nbex2=@{$coords2->{$chr}{"start"}};
	    
	    my $firstj=0;
	    
	    for(my $i=0; $i<$nbex1; $i++){
		
		my $start1=${$coords1->{$chr}{"start"}}[$i]-$margin;
		my $end1=${$coords1->{$chr}{"end"}}[$i]+$margin;
		my $strand1=${$coords1->{$chr}{"strand"}}[$i];
		my $id1=${$coords1->{$chr}{"id"}}[$i];
		
		my $j=$firstj;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}
		
		$firstj=$j;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $strand2=${$coords2->{$chr}{"strand"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];
		    
		    if(($strand1 eq $strand2 && $type eq "sense") || ($strand1 ne $strand2 && $type eq "antisense")){
			
			my $M=max($start1,$start2);
			my $m=min($end1,$end2);
			
			if($M<=$m){
			    my $fr1=($m-$M+1)/($end1-$start1+1);
			    my $fr2=($m-$M+1)/($end2-$start2+1);
			    
			    if(exists $overlap->{$id1}){
				$overlap->{$id1}{$id2}={"start"=>$M,"end"=>$m,"fractionoverlap1"=>$fr1,"fractionoverlap2"=>$fr2};
			    }
			    else{
				$overlap->{$id1}={$id2=>{"start"=>$M,"end"=>$m,"fractionoverlap1"=>$fr1,"fractionoverlap2"=>$fr2}};
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

sub addToCluster{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $indexclusters=$_[2];
    my $key=$_[3];

    if(exists $refconnected->{$key}){
                    
        my $indexcluster="NA";
	
        if(exists $indexclusters->{$key}){
            $indexcluster=$indexclusters->{$key}
        }       
              
        if($indexcluster eq "NA"){
            my $nbclusters=keys %{$refclusters};
            $indexcluster=$nbclusters+1;
            $refclusters->{$indexcluster}=[$key];
            $indexclusters->{$key}=$indexcluster;
        }
                
        foreach my $connection (keys %{$refconnected->{$key}}){
            if(!(exists $indexclusters->{$connection})){
                push(@{$refclusters->{$indexcluster}},$connection);
                $indexclusters->{$connection}=$indexcluster;
                addToCluster($refconnected,$refclusters,$indexclusters,$connection);
            }
        }

        delete $refconnected->{$key};
    }
}

###########################################################################

sub extractClusters{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    
    my %indexclusters;
    
    my $nbconnected=keys %{$refconnected};
    
    my $round=0;
    
    my %alreadyprinted;
    
    while($nbconnected>0){
	
        foreach my $key (keys %{$refconnected}){
            addToCluster($refconnected,$refclusters,\%indexclusters,$key);
            $nbconnected=keys %{$refconnected};
        }
        
        $round++;
        
        $nbconnected=keys %{$refconnected};
    }
}

##############################################################

sub readGeneInfo{
    my $pathInfo=$_[0];
    my $info=$_[1];
    
    open(my $input,$pathInfo);
    
    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);
    
    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $id=$s[$header{"stable_id"}];
	my $thisbio=$s[$header{"biotype"}];
	my $descr=$s[$header{"description"}];

	$info->{$id}={"biotype"=>$thisbio,"description"=>$descr};

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

	my $geneid=$s[$header{"GeneID"}];

	$id=$geneid.",".$id;
	
	$projectedexons->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand, "gene"=>$geneid};
		

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub constructOverlapMap{
    my $exonoverlap=$_[0];
    my $exoninfo1=$_[1];
    my $exoninfo2=$_[2];
    my $genemap=$_[3];

    my %allposmap;

    foreach my $exon1 (keys %{$exonoverlap}){
	my $gene1="NA";
	    
	if(exists $exoninfo1->{$exon1}){
	    $gene1=$exoninfo1->{$exon1}{"gene"};
	}
	else{
	    print "Weird! cannot find exon info for ".$exon1."\n";
	    exit(1);
	}

	if(!(exists $allposmap{$gene1})){
	    $allposmap{$gene1}={};
	}
	
	foreach my $exon2 (keys %{$exonoverlap->{$exon1}}){
	    my $gene2="NA";
	    
	    if(exists $exoninfo2->{$exon2}){
		$gene2=$exoninfo2->{$exon2}{"gene"};
	    }
	    else{
		print "Weird! cannot find exon info for ".$exon2."\n";
		exit(1);
	    }

	    my $startov=$exonoverlap->{$exon1}{$exon2}{"start"};
	    my $endov=$exonoverlap->{$exon1}{$exon2}{"end"};

	    for(my $i=$startov; $i<=$endov; $i++){
		if(exists $allposmap{$gene1}{$i}){
		    $allposmap{$gene1}{$i}{$gene2}=1;
		}
		else{
		    $allposmap{$gene1}{$i}={$gene2=>1};
		}
	    }
	}
    }

    foreach my $gene1 (keys %allposmap){
    	$genemap->{$gene1}={};
	
	foreach my $pos (keys %{$allposmap{$gene1}}){
	    my @ovgenes=keys %{$allposmap{$gene1}{$pos}};
	    
	    foreach my $gene2 (@ovgenes){
		if(exists $genemap->{$gene1}{$gene2}){
		    $genemap->{$gene1}{$gene2}++;
		}
		else{
		    $genemap->{$gene1}{$gene2}=1;
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
    print "This script extracts overlaps with projections.\n";
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

$parameters{"species1"}="NA";
$parameters{"species2"}="NA";
$parameters{"pathExonBlocks1"}="NA";
$parameters{"pathExonBlocks2"}="NA";
$parameters{"pathProjectedExons12"}="NA";
$parameters{"pathProjectedExons21"}="NA";
$parameters{"pathProjectionMap12"}="NA";
$parameters{"pathProjectionMap21"}="NA";
$parameters{"pathGeneClusters"}="NA";

my @defaultpars=("species1", "species2", "pathExonBlocks1", "pathExonBlocks2", "pathProjectedExons12",  "pathProjectedExons21", "pathProjectionMap12", "pathProjectionMap21", "pathGeneClusters");

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

my $sp1=$parameters{"species1"};
my $sp2=$parameters{"species2"};

print "First species: ".$sp1.".\n";
print "Second species: ".$sp2.".\n";

##############################################################

print "Reading exon coordinates...\n";

my %exons1;
my %genes1;
readExonBlocks($parameters{"pathExonBlocks1"}, \%exons1, \%genes1);

my $nbex1=keys %exons1;

print "Found ".$nbex1." exons for first species.\n";


my %exons2;
my %genes2;
readExonBlocks($parameters{"pathExonBlocks2"}, \%exons2, \%genes2);

my $nbex2=keys %exons2;

print "Found ".$nbex2." exons for second species.\n";

print "Done.\n";

##############################################################

print "Extracting gene coordinates...\n";

my %genecoords1;
extractGeneCoords(\%genes1, \%exons1, \%genecoords1);

my %genecoords2;
extractGeneCoords(\%genes2, \%exons2, \%genecoords2);

print "Done.\n";

##############################################################

print "Constructing exon blocks for each gene...\n";

my %exonblocks1;
makeExonBlocks(\%genes1, \%exons1, 0, \%exonblocks1);  

my %exonblocks2;
makeExonBlocks(\%genes2, \%exons2, 0, \%exonblocks2);  

print "Done.\n";

##############################################################

print "Ordering exons...\n";

my %orderedexons1;

orderExons(\%exons1, \%orderedexons1);

my %orderedexons2;

orderExons(\%exons2, \%orderedexons2);

print "Done.\n";

##############################################################

print "Reading projected exons...\n";

my %projected12; ## coordinates in species 2 
readProjectedExons($parameters{"pathProjectedExons12"},  \%projected12);

my $nbproj12=keys %projected12;
print "Found ".$nbproj12." projected exons from first species to second species.\n";

my %projected21; ## coordinates in species 1
readProjectedExons($parameters{"pathProjectedExons21"}, \%projected21);

my $nbproj21=keys %projected21;
print "Found ".$nbproj21." projected exons from second species to first species.\n";

print "Done.\n";

##############################################################

print "Ordering exons...\n";

my %orderedprojected12;
orderExons(\%projected12, \%orderedprojected12);

my %orderedprojected21;
orderExons(\%projected21, \%orderedprojected21);

print "Done.\n";

##############################################################

print "Computing exon overlap...\n";

my %overlap21;

extractOverlap(\%orderedprojected21, \%orderedexons1,  0, "sense", \%overlap21);

my $nbov1=keys %overlap21;

print "Found ".$nbov1." overlapping projected exons from first species to second species.\n";

my %overlap12;

extractOverlap(\%orderedprojected12, \%orderedexons2,  0, "sense", \%overlap12);


my $nbov2=keys %overlap12;

print "Found ".$nbov2." overlapping projected exons from second species to first species.\n";

##############################################################

print "Constructing exon overlap map...\n";

print "First species projected on second species\n";

my %geneoverlap12;
constructOverlapMap(\%overlap12, \%exons1, \%exons2, \%geneoverlap12);

print "Second species projected on first species\n";

my %geneoverlap21;
constructOverlapMap(\%overlap21, \%exons2, \%exons1, \%geneoverlap21);

  
##############################################################

my %connectedgenes;

print "Writing output for overlap map...\n";

open(my $output, ">".$parameters{"pathProjectionMap12"});

print $output "ID.".$sp1."\tTotalLength\tOverlap.".$sp2."\n";

foreach my $gene (keys %exonblocks1){
    my $clustid1=$sp1.":".$gene;
    
    if(!exists $connectedgenes{$clustid1}){
	$connectedgenes{$clustid1}={$clustid1=>1};
    }
   
    my $totlen=0;
    my $nbexons=@{$exonblocks1{$gene}{"start"}};

    for(my $i=0; $i<$nbexons; $i++){
	my $start=${$exonblocks1{$gene}{"start"}}[$i];
	my $end=${$exonblocks1{$gene}{"end"}}[$i];

	$totlen+=($end-$start+1);
    }

    my $line=$gene."\t".$totlen."\t";
    
    my $overlap="NA;";
   
    if(exists $geneoverlap12{$gene}){
	$overlap="";
	
	foreach my $otherg (keys %{$geneoverlap12{$gene}}){
	    my $lenov=$geneoverlap12{$gene}{$otherg};
	    $overlap.=$otherg.":".$lenov.";";
	    
	    my $clustid2=$sp2.":".$otherg;

	    $connectedgenes{$clustid1}{$clustid2}=1;
	    
	    if(!exists $connectedgenes{$clustid2}){
		$connectedgenes{$clustid2}={$clustid2=>1, $clustid1=>1};
	    }
	    else{
		$connectedgenes{$clustid2}{$clustid1}=1;
	    }    
	}
    }

    chop $overlap;
    
    $line.=$overlap."\n";

    print $output $line;
}


## second species

open(my $output, ">".$parameters{"pathProjectionMap21"});

print $output "ID.".$sp2."\tTotalLength\tOverlap.".$sp1."\n";

foreach my $gene (keys %exonblocks2){
    my $clustid1=$sp2.":".$gene;
    
    if(!exists $connectedgenes{$clustid1}){
	$connectedgenes{$clustid1}={$clustid1=>1};
    }
    
    my $totlen=0;
    my $nbexons=@{$exonblocks2{$gene}{"start"}};
    
    for(my $i=0; $i<$nbexons; $i++){
	my $start=${$exonblocks2{$gene}{"start"}}[$i];
	my $end=${$exonblocks2{$gene}{"end"}}[$i];

	$totlen+=($end-$start+1);
    }


    my $line=$gene."\t".$totlen."\t";

    my $overlap="NA;";
   
    if(exists $geneoverlap21{$gene}){
	$overlap="";
	
	foreach my $otherg (keys %{$geneoverlap21{$gene}}){
	    my $lenov=$geneoverlap21{$gene}{$otherg};

	    $overlap.=$otherg.":".$lenov.";";

	    my $clustid2=$sp1.":".$otherg;
	    
	    $connectedgenes{$clustid1}{$clustid2}=1;
	    
	    if(!exists $connectedgenes{$clustid2}){
		$connectedgenes{$clustid2}={$clustid2=>1, $clustid1=>1};
	    } else{
		$connectedgenes{$clustid2}{$clustid1}=1;
	    }
	}
    }

    chop $overlap;
    
    $line.=$overlap."\n";

    print $output $line;
}

print "Done.\n";

##############################################################

print "Extracting gene clusters...\n";

my %geneclusters;
extractClusters(\%connectedgenes, \%geneclusters);

print "Done.\n";

##############################################################

print "Writing output for clusters...\n";

open(my $output, ">".$parameters{"pathGeneClusters"});

print $output "Genes.".$sp1."\tGenes.".$sp2."\n";

foreach my $indexcluster (keys %geneclusters){
    my %genes1;
    my %genes2;

    foreach my $gene (@{$geneclusters{$indexcluster}}){
	my @s=split(":", $gene);

	if($s[0] eq $sp1){
	    $genes1{$s[1]}=1;
	} else{
	    if($s[0] eq $sp2){
		$genes2{$s[1]}=1;
	    }  else{
		print "Weird! cannot recognize species: ".$s[0]."\n";
		exit(1);
	    }
	}
    }

    my $id1="NA";
    my $id2="NA";

    my $nb1=keys %genes1;
    
    if($nb1>0){
	$id1=join(";", keys %genes1);
    }

    my $nb2=keys %genes2;
    
    if($nb2>0){
	$id2=join(";", keys %genes2);
    }

    print $output $id1."\t".$id2."\n";
}

close($output);

print "Done.\n";

##############################################################
