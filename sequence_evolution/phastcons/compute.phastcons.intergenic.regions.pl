use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readPhastCons{
    my $pathin=$_[0];
    my $refphast=$_[1];
    
    my @s=split("\\.",$pathin);
    my $nbs=@s;
    my $ext=$s[$nbs-1];

    my $input;
    
    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    }
    else{
	open($input,$pathin);
    }
    
    my $line=<$input>;
    my $currentpos="NA";
    
    while($line){
	my $first=substr $line,0,1;
	
	if($first eq "f"){
	    my @s=split(" ",$line);
	    my $start=$s[2];
	    my @t=split("=",$start);
	    $currentpos=$t[1]+0;
	    $line=<$input>;
	    next;
	}
	else{
	    chomp $line;
	    my $val=$line+0.0;

	    $refphast->{$currentpos}=$val;
	    
	    $currentpos++;
	    $line=<$input>;
	}
    }
    
    close($input);
}

##############################################################

sub computePhastConsScore{
    my $coords=$_[0];
    my $phast=$_[1];
    my $maskedcoords=$_[2];
    my $chr=$_[3];

    my $nbmotifs=@{$coords->{$chr}{"start"}};

    for(my $i=0;$i<$nbmotifs;$i++){
	my $start=${$coords->{$chr}{"start"}}[$i];
	my $end=${$coords->{$chr}{"end"}}[$i];

	my $sumscore=0;
	my $nbpos=0;
	my $nbunmasked=0;

	for(my $j=$start;$j<=$end;$j++){
	    if(!(exists $maskedcoords->{$j})){
		$nbunmasked++;
		if(exists $phast->{$j}){
		    $nbpos++;
		    $sumscore+=$phast->{$j};
		}
	    }
	}
	
	if($nbpos>0){
	    ${$coords->{$chr}{"score"}}[$i]=($sumscore+0.0)/($nbpos+0.0);
	}
	
	${$coords->{$chr}{"coveredbases"}}[$i]=$nbpos+0.0;
	${$coords->{$chr}{"unmaskedbases"}}[$i]=$nbunmasked+0.0;
    }
    
}

##############################################################

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

sub constructHashCoordinates{
    my $refblocks=$_[0];
    my $chr=$_[1];
    my $hashcoords=$_[2];

    foreach my $gene (keys %{$refblocks}){
	my $ens=$refblocks->{$gene}{"chr"};
	my $ucsc="chr".$ens;

	if($chr eq $ucsc || $chr eq $ens){
	    my $nbex=@{$refblocks->{$gene}{"start"}};
	    
	    for(my $i=0; $i<$nbex; $i++){
		my $start=${$refblocks->{$gene}{"start"}}[$i];
		my $end=${$refblocks->{$gene}{"end"}}[$i];

		for(my $pos=$start; $pos<=$end; $pos++){
		    $hashcoords->{$pos}=1;
		}
	    }
	}
    }
}

##########################################################################

sub readCoords{
    my $pathin=$_[0];
    my $refcoords=$_[1];
    
    open(my $input, $pathin);
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $chr=$s[0];

	my $prefix=substr $chr,0,3;

	if($prefix ne "chr"){
	    $chr="chr".$chr;
	}
	
	my $start=$s[1]+0; ## 1-based
	my $end=$s[2]+0;
	my $strand="+";
	my $id=$chr.",".$start.",".$end;

	$refcoords->{$id}={"chr"=>$chr, "strand"=>$strand, "start"=>[$start], "end"=>[$end]};
	 	
	$line=<$input>;
    }
    
    close($input);
}

##########################################################################

sub orderCoords{
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
	$refordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"gene"=>[],"score"=>[],"coveredbases"=>[],"unmaskedbases"=>[]};

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
			push(@{$refordered->{$chr}{"score"}},"NA");
			push(@{$refordered->{$chr}{"coveredbases"}},0);
			push(@{$refordered->{$chr}{"unmaskedbases"}},0);
		    }
		} 
	    }
	} 
    }

}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes coverage for promoters. \n";
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

$parameters{"pathIntergenicRegionCoords"}="NA";
$parameters{"pathPhastCons"}="NA";
$parameters{"chr"}="NA";
$parameters{"pathOutput"}="NA";


my %defaultvalues;
my @defaultpars=("pathIntergenicRegionCoords",  "pathPhastCons", "chr", "pathOutput");

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


##############################################################
##############################################################

my %maskexons;
my %maskedcoords;

##############################################################

print "Reading intergenic regions...\n";
my %intergenic;
readCoords($parameters{"pathIntergenicRegionCoords"}, \%intergenic);
print "Done.\n";

print "Ordering intergenic...\n";
my %orderedintergenic;
orderCoords(\%intergenic, \%orderedintergenic);
print "Done\n";

print "Computing PhastCons score and writing output \n";

open(my $output,">".$parameters{"pathOutput"});
print $output "ID\tChr\tStart\tEnd\tStrand\tScore\tCoveredLength\tUnmaskedLength\n";

my $chr=$parameters{"chr"};
my $path=$parameters{"pathPhastCons"};

print "Chromosome ".$chr."\n";


if(-e $path){
    if(exists $orderedintergenic{$chr}){
	
	my %phastCons;
	readPhastCons($path,\%phastCons);
	
	computePhastConsScore(\%orderedintergenic,\%phastCons,\%maskedcoords, $chr);
        
	my $nborderedintergenic=@{$orderedintergenic{$chr}{"start"}};
	
	my %scoregenes;
	
	
	for(my $i=0;$i<$nborderedintergenic;$i++){
	    
	    my $gene=${$orderedintergenic{$chr}{"gene"}}[$i];
	    
	    print $output ${$orderedintergenic{$chr}{"gene"}}[$i]."\t".$chr."\t".${$orderedintergenic{$chr}{"start"}}[$i]."\t".${$orderedintergenic{$chr}{"end"}}[$i]."\t".${$orderedintergenic{$chr}{"strand"}}[$i]."\t".${$orderedintergenic{$chr}{"score"}}[$i]."\t".${$orderedintergenic{$chr}{"coveredbases"}}[$i]."\t".${$orderedintergenic{$chr}{"unmaskedbases"}}[$i]."\n";
	    
	}
    }
}

close($output);

print "Done.\n";


##############################################################
