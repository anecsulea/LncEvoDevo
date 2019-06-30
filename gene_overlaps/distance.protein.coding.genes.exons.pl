use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readGTF{
    my $pathin=$_[0];
    my $forbiddentx=$_[1];
    my $selectedgenes=$_[2];  
    my $genecoords=$_[3];
    my $exoncoords=$_[4];
    
    my $nbf=keys %{$forbiddentx};
    print "There are ".$nbf." forbidden transcripts.\n";

    my $nbsel=keys %{$selectedgenes};
    my $activegenesel="no";
    
    if($nbsel>0){
	print "There are ".$nbsel." selected genes.\n";
	$activegenesel="yes";
    }
    else{
	print "We keep all annotated genes.\n";
    }
    
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

		if(exists $forbiddentx->{$tx}){
		    print "Discarding transcript ".$tx."\n";
		}
		else{
		    if($gene eq "NA"){
			print "Weird! cannot find gene id for ".$line."\n";
			exit(1);
		    }

		    if($activegenesel eq "no" || exists $selectedgenes->{$gene}){

			my $exonid=$chr.",".$start.",".$end.",".$strand;

			$exoncoords->{$exonid}={"start"=>$start, "end"=>$end, "chr"=>$chr, "strand"=>$strand};
		    
			if(exists $genecoords->{$gene}){
			    if($start < $genecoords->{$gene}{"start"}){
				$genecoords->{$gene}{"start"}=$start;
			    }
			    
			    if($end > $genecoords->{$gene}{"end"}){
				$genecoords->{$gene}{"end"}=$end;
			    }
			    
			} else{
			    $genecoords->{$gene}={"start"=>$start, "end"=>$end, "chr"=>$chr, "strand"=>$strand};
			}
		    }
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

sub orderCoordsStart{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $exid (keys %{$coords}){
	my $chr=$coords->{$exid}{"chr"};
	my $strand=$coords->{$exid}{"strand"};
	my $start=$coords->{$exid}{"start"};
	my $end=$coords->{$exid}{"end"};

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
	$ordered->{$chr}={};

	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashpos{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		my $end=${$hashpos{$chr}{$start}{"end"}}[$i];
		my $strand=${$hashpos{$chr}{$start}{"strand"}}[$i];

		if(!(exists $ordered->{$chr}{$strand})){
		    $ordered->{$chr}{$strand}={"start"=>[],"end"=>[],"strand"=>[],"id"=>[]};
		}

		my $id=${$hashpos{$chr}{$start}{"id"}}[$i];
	
		push(@{$ordered->{$chr}{$strand}{"start"}}, $start);
		push(@{$ordered->{$chr}{$strand}{"end"}}, $end);
		push(@{$ordered->{$chr}{$strand}{"id"}}, $id);
	    }
	}
    }
}

##############################################################################

sub orderCoordsEnd{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $exid (keys %{$coords}){
	my $chr=$coords->{$exid}{"chr"};
	my $strand=$coords->{$exid}{"strand"};
	my $start=$coords->{$exid}{"start"};
	my $end=$coords->{$exid}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$end}){
		push(@{$hashpos{$chr}{$end}{"id"}},$exid);
		push(@{$hashpos{$chr}{$end}{"start"}},$start);
		push(@{$hashpos{$chr}{$end}{"strand"}},$strand);
	    }
	    else{
		$hashpos{$chr}{$end}={"id"=>[$exid],"start"=>[$start],"strand"=>[$strand]};
	    }
	}
	else{
	    $hashpos{$chr}={$end=>{"id"=>[$exid],"start"=>[$start],"strand"=>[$strand]}};
	}
    }

    foreach my $chr (keys %hashpos){
	$ordered->{$chr}={};

	my @uniqueend=keys %{$hashpos{$chr}};
	my @sortedend=sort {$a<=>$b} @uniqueend;

	foreach my $end (@sortedend){
	    my $nb=@{$hashpos{$chr}{$end}{"start"}};

	    for(my $i=0; $i<$nb; $i++){
		my $start=${$hashpos{$chr}{$end}{"start"}}[$i];
		my $strand=${$hashpos{$chr}{$end}{"strand"}}[$i];
		my $id=${$hashpos{$chr}{$end}{"id"}}[$i];

		if(!(exists $ordered->{$chr}{$strand})){
		    $ordered->{$chr}{$strand}={"start"=>[],"end"=>[],"strand"=>[],"id"=>[]};
		}
	
		push(@{$ordered->{$chr}{$strand}{"start"}}, $start);
		push(@{$ordered->{$chr}{$strand}{"end"}}, $end);
		push(@{$ordered->{$chr}{$strand}{"id"}}, $id);
	    }
	}
    }
}

###########################################################################

sub computeDistanceLeft{
    my $orderedcoords=$_[0]; ## they have to be ordered by the end coordinate
    my $geneids=$_[1];
    my $exonids=$_[2];
    my $distance=$_[3];

    foreach my $chr (keys %{$orderedcoords}){
	foreach my $strand (keys %{$orderedcoords->{$chr}}){
	    my $nb=@{$orderedcoords->{$chr}{$strand}{"start"}};

	    for(my $i=0; $i<$nb; $i++){
		my $id=${$orderedcoords->{$chr}{$strand}{"id"}}[$i];
		
		if(exists $geneids->{$id}){
		    ## we need to find the closest exon
		    my $start=${$orderedcoords->{$chr}{$strand}{"start"}}[$i];
		    my $end=${$orderedcoords->{$chr}{$strand}{"end"}}[$i];
		    
		    my $dist="Inf";
		    my $closestexon="NA";

		    my $j=$i-1;
		    
		    while($j>=0){
			my $nextid=${$orderedcoords->{$chr}{$strand}{"id"}}[$j];
			
			if(exists $exonids->{$nextid}){
			    my $nextstart=${$orderedcoords->{$chr}{$strand}{"start"}}[$j];
			    my $nextend=${$orderedcoords->{$chr}{$strand}{"end"}}[$j];
			    my $nextdist=$start-$nextend;

			    if($nextdist<0){
				$nextdist=-1; ## we homogenize all negative distances to a value of -1
			    }
			    
			    $dist=$nextdist;
			    $closestexon=$nextid;
			    
			    last; ## we stop at the first exon we see
			}
			
			$j--;
		    }
		    
		    $distance->{$id}={"distance"=>$dist, "closest"=>$closestexon};
		}
	    }
	}
    }
}

###########################################################################

sub computeDistanceRight{
    my $orderedcoords=$_[0]; ## they have to be ordered by the start coordinate
    my $geneids=$_[1];
    my $exonids=$_[2];
    my $distance=$_[3];

    foreach my $chr (keys %{$orderedcoords}){
	foreach my $strand (keys %{$orderedcoords->{$chr}}){
	    my $nb=@{$orderedcoords->{$chr}{$strand}{"start"}};

	    for(my $i=0; $i<$nb; $i++){
		my $id=${$orderedcoords->{$chr}{$strand}{"id"}}[$i];
		
		if(exists $geneids->{$id}){
		    ## we need to find the closest exon
		    my $start=${$orderedcoords->{$chr}{$strand}{"start"}}[$i];
		    my $end=${$orderedcoords->{$chr}{$strand}{"end"}}[$i];
		    
		    my $dist="Inf";
		    my $closestexon="NA";

		    my $j=$i+1;
		    
		    while($j<$nb){
			my $nextid=${$orderedcoords->{$chr}{$strand}{"id"}}[$j];
			
			if(exists $exonids->{$nextid}){
			    my $nextstart=${$orderedcoords->{$chr}{$strand}{"start"}}[$j];
			    my $nextend=${$orderedcoords->{$chr}{$strand}{"end"}}[$j];
			    my $nextdist=$nextstart-$end;

			    if($nextdist<0){
				$nextdist=-1; ## we homogenize all negative distances to a value of -1
			    }
			    
			    $dist=$nextdist;
			    $closestexon=$nextid;
			    
			    last; ## we stop at the first exon we see
			}
			
			$j++;
		    }
		    
		    $distance->{$id}={"distance"=>$dist, "closest"=>$closestexon};
		}
	    }
	}
    }
}


##############################################################################

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

###########################################################################

sub readSynonyms{
    my $pathin=$_[0];
    my $syn=$_[1];

    open(my $input, $pathin);

    my $line=<$input>; ## header
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my @oldids=split(",", $s[1]);
	my $newid=$s[2];

	if(!(exists $syn->{$newid})){
	   $syn->{$newid}={};
	}
	
	foreach my $oldid (@oldids){
	    $syn->{$newid}{$oldid}=1;
	}
	
	$line=<$input>;
    }

    close($input);
}

###########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes the distance between newly annotated loci and protein-coding gene exons. \n";
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

$parameters{"pathGTF1"}="NA";
$parameters{"pathGTF2"}="NA";
$parameters{"forbiddenTranscripts"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"pathSynonyms"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathGTF1", "pathGTF2", "forbiddenTranscripts", "pathGeneInfo", "pathSynonyms", "pathOutput");

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

my %exoncoords1;
my %genecoords1;

my %forbiddentx1;
my %selectedgenes1;

print "There are no forbidden transcripts or selected genes for first annotation set.\n";

readGTF($parameters{"pathGTF1"}, \%forbiddentx1, \%selectedgenes1, \%genecoords1, \%exoncoords1);
  
my $nbg1=keys %genecoords1;
my $nbex1=keys %exoncoords1;

print "Found ".$nbg1." genes and ".$nbex1." exons in the first annotation set.\n";

print "Done.\n";

#####################################################################

print "Reading gene info for second annotation set...\n";

my %geneinfo;
readGeneInfo($parameters{"pathGeneInfo"}, \%geneinfo);

my $nbg=keys %geneinfo;

print "Found ".$nbg." genes with info in the second annotation set.\n";

print "Done.\n";

print "Reading gene synoynms...\n";

my %synonyms;

readSynonyms($parameters{"pathSynonyms"}, \%synonyms);

my $nbsyn=keys %synonyms;

print "Done.\n";


print "Extracting protein-coding genes...\n";

my %pcgenes;

foreach my $gene (keys %geneinfo){
    my $bio=$geneinfo{$gene}{"biotype"};
    
    if($bio eq "protein_coding"){
	$pcgenes{$gene}=1;
    }
}

foreach my $gene (keys %synonyms){
    foreach my $oldid (keys %{$synonyms{$gene}}){
	if(exists $geneinfo{$oldid} && $geneinfo{$oldid}{"biotype"} eq "protein_coding"){
	    $pcgenes{$gene}=1;
	}
    }
}

my $nbpc=keys %pcgenes;

print "Found ".$nbpc." protein-coding genes.\n";

print "Done.\n";

#####################################################################

print "Reading coordinates for second annotation set...\n";

my %forbiddentx2;
my @s=split(",",$parameters{"forbiddenTranscripts"});

foreach my $f (@s){
    $forbiddentx2{$f}=1;
}

my $nbf=keys %forbiddentx2;

print "There are ".$nbf." forbidden transcripts.\n";

my %genecoords2;
my %exoncoords2;

readGTF($parameters{"pathGTF2"}, \%forbiddentx2, \%pcgenes, \%genecoords2, \%exoncoords2);

my $nbex2=keys %exoncoords2;
my $nbg2=keys %genecoords2;

print "Found ".$nbg2." genes and ".$nbex2." exons in second annotation set (protein-coding only)\n";
 
print "Done.\n";

#####################################################################

print "Combining gene and exon coordinates...\n";

my %allcoords;

foreach my $gene (keys %genecoords1){
    $allcoords{$gene}={};
    foreach my $key (keys %{$genecoords1{$gene}}){
	$allcoords{$gene}{$key}=$genecoords1{$gene}{$key};
    }
}

foreach my $exon (keys %exoncoords2){
    $allcoords{$exon}={};
    foreach my $key (keys %{$exoncoords2{$exon}}){
	$allcoords{$exon}{$key}=$exoncoords2{$exon}{$key};
    }
}

print "Done.\n";

#####################################################################

print "Ordering coordinates...\n";

my %orderedcoordsstart;
orderCoordsStart(\%allcoords, \%orderedcoordsstart);

my %orderedcoordsend;
orderCoordsEnd(\%allcoords, \%orderedcoordsend);

print "Done.\n";

#####################################################################

print "Computing distances...\n";

my %distanceleft;
computeDistanceLeft(\%orderedcoordsend, \%genecoords1, \%exoncoords2, \%distanceleft);

my $nbdl=keys %distanceleft;

print "computed distance left for ".$nbdl." genes.\n";
 
my %distanceright;
computeDistanceRight(\%orderedcoordsstart, \%genecoords1, \%exoncoords2, \%distanceright);

my $nbdr=keys %distanceright;

print "computed distance right for ".$nbdr." genes.\n";
 
 print "Done.\n";

#####################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tChr\tStart\tEnd\tStrand\tDistanceLeft\tClosestIDLeft\tDistanceRight\tClosestIDRight\n";

foreach my $gene (keys %genecoords1){
    my $chr=$genecoords1{$gene}{"chr"};
    my $start=$genecoords1{$gene}{"start"};
    my $end=$genecoords1{$gene}{"end"};
    my $strand=$genecoords1{$gene}{"strand"};

    if(exists $distanceleft{$gene}){
	if(exists $distanceright{$gene}){
	    my $distleft=$distanceleft{$gene}{"distance"};
	    my $closestleft=$distanceleft{$gene}{"closest"};

	    my $distright=$distanceright{$gene}{"distance"};
	    my $closestright=$distanceright{$gene}{"closest"};

	    print $output $gene."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$distleft."\t".$closestleft."\t".$distright."\t".$closestright."\n";
	}	
	else{
	    print "Weird! cannot find distance right for ".$gene."\n";
	    exit(1);
	}
    }
    else{
	print "Weird! cannot find distance left for ".$gene."\n";
	exit(1);
    }
}

close($output);

print "Done.\n";

#####################################################################

