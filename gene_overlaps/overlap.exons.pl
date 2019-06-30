use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $refblocks=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $gene=$s[0];
	my $idexon=$s[1];
	my $chr=$s[2];
	my $start=$s[3]+0;
	my $end=$s[4]+0;
	my $strand=$s[5];

	if($strand eq "."){
	    $strand="NA";
	}

	if(exists $refblocks->{$gene}){
	    push(@{$refblocks->{$gene}{"start"}},$start);
	    push(@{$refblocks->{$gene}{"end"}},$end);
	}
	else{
	    $refblocks->{$gene}={"chr"=>$chr,"strand"=>$strand,"start"=>[$start],"end"=>[$end]};
	}
	    
	$line=<$input>;
    }
    close($input);
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


##############################################################

sub overlapBlocks{
    my $refblocks1=$_[0]; ## ordered exon blocks
    my $refblocks2=$_[1]; ## ordered exon blocks
    my $type=$_[2];
    my $refoverlap=$_[3]; ## overlap with blocks
   
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
		    
		    my $strand2=${$refblocks2->{$chr}{"strand"}}[$j];

		    if(($type eq "sense" && $strand1 eq $strand2) || ($type eq "antisense" && $strand1 ne $strand2) || ($type eq "any")){
			my $lenov=$m-$M+1;
			
			if($lenov>=1){
			    if(exists $refoverlap->{$gene1}){
				push(@{$refoverlap->{$gene1}{"start"}},$M);
				push(@{$refoverlap->{$gene1}{"end"}},$m);
			    }
			    else{
				$refoverlap->{$gene1}={"start"=>[$M],"end"=>[$m]};
			    }
			}
		    }
		    
		    $j++;
		}
	    }
	}
    }
}

##############################################################

sub makeBlocks{
    my $intervals=$_[0];
    my $blocks=$_[1];

    my %hashpos;

    my $nb=@{$intervals->{"start"}};
    $blocks->{"start"}=[];
    $blocks->{"end"}=[];

    for(my $i=0; $i<$nb; $i++){
	my $start=${$intervals->{"start"}}[$i];
	my $end=${$intervals->{"end"}}[$i];

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
    my $totlen=0;

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
	    push(@{$blocks->{"start"}}, $currentstart);
	    push(@{$blocks->{"end"}}, $currentend);
	    $totlen+=($currentend-$currentstart+1);

	    $currentstart=$thisstart;
	    $currentend=$thisend;
	}
    }
    
    push(@{$blocks->{"start"}}, $currentstart);
    push(@{$blocks->{"end"}}, $currentend);
    $totlen+=($currentend-$currentstart+1);

    $blocks->{"totallength"}=$totlen;
    
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

############################################################################

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

################################################################################
################################################################################

my %parameters;

$parameters{"pathExonBlocks1"}="NA";
$parameters{"pathExonBlocks2"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"biotypes"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonBlocks1", "pathExonBlocks2", "pathGeneInfo", "biotypes", "pathOutput");

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

print "Reading exon blocks for first annotation set...\n";

my %exonblocks1;

readExonBlocks($parameters{"pathExonBlocks1"}, \%exonblocks1);

my $nbg=keys %exonblocks1;

print "Found ".$nbg." genes in first annotation set.\n";

print "Done.\n";


print "Reading exon blocks for second annotation set...\n";

my %exonblocks2;

readExonBlocks($parameters{"pathExonBlocks2"}, \%exonblocks2);

my $nbg=keys %exonblocks2;

print "Found ".$nbg." genes in second annotation set.\n";

print "Done.\n";

#####################################################################

print "Reading gene info...\n";

my %geneinfo;

readGeneInfo($parameters{"pathGeneInfo"}, \%geneinfo);

print "Done.\n";

#####################################################################

print "Filtering genes from the second annotation set...\n";

my @biotypes=split(",", $parameters{"biotypes"});
my %biotypes;

foreach my $b (@biotypes){
    $biotypes{$b}=1;
}
my $nb=keys %biotypes;

print "There are ".$nb." accepted biotypes: ".join("; ",keys %biotypes)."\n";

my @allgenes=keys %exonblocks2;
my $nbkept=0;

foreach my $gene (@allgenes){
    if(exists $geneinfo{$gene} && exists $biotypes{$geneinfo{$gene}{"biotype"}}){
	$nbkept++;
    } else{
	delete $exonblocks2{$gene};
    }
}

my $nbrealkept=keys %exonblocks2;

if($nbrealkept != $nbkept){
    print "Weird! saw different numbers of kept genes: ".$nbkept." ".$nbrealkept."\n";
    exit(1);
} 

print "Done.\n";

#####################################################################

print "Ordering exon blocks...\n";

my %orderedblocks1;
orderExonBlocks(\%exonblocks1, \%orderedblocks1);

my %orderedblocks2;
orderExonBlocks(\%exonblocks2, \%orderedblocks2);

print "Done.\n";

#####################################################################

print "Extracting overlap...\n";

my %overlapsense;
my %overlapantisense;
my %overlapboth;

overlapBlocks(\%orderedblocks1, \%orderedblocks2, "sense", \%overlapsense);

overlapBlocks(\%orderedblocks1, \%orderedblocks2, "antisense", \%overlapantisense);

overlapBlocks(\%orderedblocks1, \%orderedblocks2, "any", \%overlapboth);

my $nbovsense=keys %overlapsense;
my $nbovantisense=keys %overlapantisense;
my $nbovboth=keys %overlapboth;

print "Found ".$nbovsense." genes with overlap on the sense strand.\n";
print "Found ".$nbovantisense." genes with overlap on the antisense strand.\n";
print "Found ".$nbovboth." genes with overlap on both strands.\n";

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tTotalExonicLength\tLengthOverlapSense\tLengthOverlapAntisense\tLengthOverlapBothStrands\n";

foreach my $gene (keys %exonblocks1){
    my $totlen=0;
    
    my $nbex=@{$exonblocks1{$gene}{"start"}};
    
    for(my $i=0; $i<$nbex; $i++){
	$totlen+=(${$exonblocks1{$gene}{"end"}}[$i] - ${$exonblocks1{$gene}{"start"}}[$i]+1);
    }

    my $lensense=0;
    my $lenantisense=0;
    my $lenboth=0;

    if(exists $overlapsense{$gene}){
	my %blocks;
	makeBlocks($overlapsense{$gene}, \%blocks);
	my $nbb=@{$blocks{"start"}};

	for(my $i=0; $i<$nbb; $i++){
	    $lensense+=(${$blocks{"end"}}[$i] - ${$blocks{"start"}}[$i]+1);
	}
    }

    if(exists $overlapantisense{$gene}){
	my %blocks;
	makeBlocks($overlapantisense{$gene}, \%blocks);
	my $nbb=@{$blocks{"start"}};

	for(my $i=0; $i<$nbb; $i++){
	    $lenantisense+=(${$blocks{"end"}}[$i] - ${$blocks{"start"}}[$i]+1);
	}
    }

    if(exists $overlapboth{$gene}){
	my %blocks;
	makeBlocks($overlapboth{$gene}, \%blocks);
	my $nbb=@{$blocks{"start"}};
	
	for(my $i=0; $i<$nbb; $i++){
	    $lenboth+=(${$blocks{"end"}}[$i] - ${$blocks{"start"}}[$i]+1);
	}
    }

    print $output $gene."\t".$totlen."\t".$lensense."\t".$lenantisense."\t".$lenboth."\n";
}

close($output);

print "Done.\n";

##############################################################
