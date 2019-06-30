use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readGTF{
    my $pathin=$_[0];
    my $forbidden=$_[1];
    my $genecoords=$_[2];
    
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

		if(exists $forbidden->{$tx}){
		    print "Discarding transcript ".$tx."\n";
		}
		else{
		    if($gene eq "NA"){
			print "Weird! cannot find gene id for ".$line."\n";
			exit(1);
		    }
		    
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

sub orderCoords{
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

###################################################################################################

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
		    
		    if(($strand1 eq $strand2 && $type eq "sense") || ($strand1 ne $strand2 && $type eq "antisense") || ($type eq "any")){
			
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

####################################################################################################

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

###############################################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes overlaps between two sets of genes. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

###############################################################################################################
###############################################################################################################

my %parameters;

$parameters{"pathGTF1"}="NA";
$parameters{"pathGTF2"}="NA";
$parameters{"forbiddenTranscripts"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"biotypes"}="NA";
$parameters{"windowSizes"}="NA";
$parameters{"overlapSense"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathGTF1", "pathGTF2", "forbiddenTranscripts", "pathGeneInfo", "biotypes", "windowSizes", "overlapSense", "pathOutput");

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

print "Reading gene coordinates...\n";

my %forbidden;
my @s=split(",", $parameters{"forbiddenTranscripts"});
foreach my $tx (@s){
    $forbidden{$tx}=1;
}

my $nbf=keys %forbidden;

if($nbf>0){
    print "There are ".$nbf." forbidden transcripts: ".join("; ",keys %forbidden)."\n";
}

my %genecoords1;
readGTF($parameters{"pathGTF1"}, \%forbidden, \%genecoords1);

my %genecoords2;
readGTF($parameters{"pathGTF2"}, \%forbidden, \%genecoords2);

my $nb1=keys %genecoords1;
my $nb2=keys %genecoords2;

print "Found ".$nb1." genes in first annotation set.\n";
print "Found ".$nb2." genes in second annotation set.\n";

print "Done.\n";

#####################################################################

print "Ordering coordinates...\n";

my %orderedgenes1;
orderCoords(\%genecoords1, \%orderedgenes1);

my %orderedgenes2;
orderCoords(\%genecoords2, \%orderedgenes2);

print "Done.\n";

#####################################################################

print "Reading gene biotypes...\n";

my %geneinfo;

readGeneInfo($parameters{"pathGeneInfo"}, \%geneinfo);
my $nbg=keys %geneinfo;
print "Found ".$nbg." genes with biotype info.\n";

print "Done.\n";

print "Extracting accepted biotypes...\n";

my @b=split(",", $parameters{"biotypes"});
my $allbio=$parameters{"biotypes"}; ## if we take all or NA biotypes
if($allbio eq "any" || $allbio eq "NA" || $allbio eq ""){
    print "We keep all biotypes.\n";
}

my %biotypes;

foreach my $bb (@b){
    $biotypes{$bb}=1;
}
my $nb=keys %biotypes;

print "There are ".$nb." accepted biotypes: ".join(";", keys %biotypes)."\n";

print "Done.\n";

#####################################################################

print "Reading window sizes...\n";

my @w=split(",", $parameters{"windowSizes"});
my %winsizes;

foreach my $ww (@w){
    my $www=$ww+0;
    $winsizes{$www}=1;
}

my @wins=keys %winsizes;
my @sortedwins=sort {$a<=>$b} @wins;

my $nbw=keys %winsizes;

print "There are ".$nbw." window sizes to be tested: ".join(";", keys %winsizes)."\n";

print "Done.\n";

#####################################################################

print "Computing overlap between genes...\n";

my $ovsense=$parameters{"overlapSense"};

print "Overlap sense: ".$ovsense."\n";

my %overlap;

foreach my $win (@sortedwins){
    $overlap{$win}={};

    extractOverlap(\%orderedgenes1, \%orderedgenes2, $win, $ovsense, $overlap{$win}); 

    my $nbgov=keys %{$overlap{$win}};

    print "Found ".$nbgov." genes with overlap.\n";
}

print "Done.\n";

#####################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

my $line="GeneID\tChr\tStart\tEnd\tStrand";

foreach my $win (@sortedwins){
    my $kb=$win/1000;
    $line.="\tOverlap".$kb."kb";
}

print $output $line."\n";

foreach my $gene (keys %genecoords1){
    my $chr=$genecoords1{$gene}{"chr"};
    my $start=$genecoords1{$gene}{"start"};
    my $end=$genecoords1{$gene}{"end"};
    my $strand=$genecoords1{$gene}{"strand"};

    my $line=$gene."\t".$chr."\t".$start."\t".$end."\t".$strand;

    foreach my $win (@sortedwins){
	my %keptov;

	if(exists $overlap{$win}{$gene}){
	    foreach my $gene2 (keys %{$overlap{$win}{$gene}}){
		if($gene2 ne $gene){
		    if($allbio eq "any" || $allbio eq "NA" || $allbio eq ""){
			$keptov{$gene2}=1;
		    }
		    else{
			if(exists $geneinfo{$gene2}){
			    my $thisb=$geneinfo{$gene2}{"biotype"};
			    
			    if(exists $biotypes{$thisb}){
				$keptov{$gene2}=1;
			    }

			} else{
			    print "Weird!! cannot find biotype for ".$gene2.", something went wrong.\n";
			    exit(1);
			}
		    }
		}
	    }
	}

	my $nbkept=keys %keptov;
	
	if($nbkept>0){
	    $line.="\t".join(";",keys %keptov);
	} else{
	    $line.="\tNA";
	}
    }
    
    print $output $line."\n";
}

close($output);

print "Done.\n";

#####################################################################
