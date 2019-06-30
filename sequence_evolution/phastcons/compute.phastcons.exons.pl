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
    my $refphylop=$_[1];
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

		if(exists $refphylop->{$j}){
		    $nbpos++;
		    $sumscore+=$refphylop->{$j};
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
	if($chr eq $refblocks->{$gene}{"chr"}){
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
    print "This script computes coverage for exons. \n";
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
$parameters{"pathMaskExonBlocks"}="NA";
$parameters{"pathPhastCons"}="NA";
$parameters{"chr"}="NA";
$parameters{"pathOutputExons"}="NA";
$parameters{"pathOutputGenes"}="NA";


my %defaultvalues;
my @defaultpars=("pathExonBlocks", "pathMaskExonBlocks","pathPhastCons","chr","pathOutputExons","pathOutputGenes");

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

my $mask="no";

if(-e $parameters{"pathMaskExonBlocks"}){
    $mask="yes";

    print "Reading masked coordinates...\n";
    readExonBlocks($parameters{"pathMaskExonBlocks"}, \%maskexons);
    print "Done.\n";
}

##############################################################

print "Reading exons...\n";
my %exons;
readExonBlocks($parameters{"pathExonBlocks"},\%exons);
print "Done.\n";

print "Ordering exons...\n";
my %orderedexons;
orderExonBlocks(\%exons,\%orderedexons);
print "Done\n";

print "Computing PhastCons score and writing output \n";

open(my $outputex,">".$parameters{"pathOutputExons"});
print $outputex "Gene\tChr\tStart\tEnd\tStrand\tScore\tCoveredLength\tUnmaskedLength\n";


open(my $outputgene,">".$parameters{"pathOutputGenes"});
print $outputgene "Gene\tChr\tStart\tEnd\tStrand\tScore\tCoveredLength\tUnmaskedLength\n";


my $chr=$parameters{"chr"};
my $path=$parameters{"pathPhastCons"};

print "Chromosome ".$chr."\n";

if($mask eq "yes"){
    print "Masking exon coordinates.\n";
    constructHashCoordinates(\%maskexons, $chr, \%maskedcoords);
}
else{
    print "Not masking anything.\n";
}

if(-e $path){
    if(exists $orderedexons{$chr}){
	
	my %phastCons;
	readPhastCons($path,\%phastCons);
	
	computePhastConsScore(\%orderedexons,\%phastCons, \%maskedcoords, $chr);
        
	my $nborderedexons=@{$orderedexons{$chr}{"start"}};
	
	my %scoregenes;
	
	
	for(my $i=0;$i<$nborderedexons;$i++){
	    
	    my $gene=${$orderedexons{$chr}{"gene"}}[$i];
	    
	    print $outputex ${$orderedexons{$chr}{"gene"}}[$i]."\t".$chr."\t".${$orderedexons{$chr}{"start"}}[$i]."\t".${$orderedexons{$chr}{"end"}}[$i]."\t".${$orderedexons{$chr}{"strand"}}[$i]."\t".${$orderedexons{$chr}{"score"}}[$i]."\t".${$orderedexons{$chr}{"coveredbases"}}[$i]."\t".${$orderedexons{$chr}{"unmaskedbases"}}[$i]."\n";
	    
	    if(exists $scoregenes{$gene}){
		$scoregenes{$gene}{"sumscore"}+=${$orderedexons{$chr}{"score"}}[$i]*${$orderedexons{$chr}{"coveredbases"}}[$i];
		$scoregenes{$gene}{"exoniclength"}+=${$orderedexons{$chr}{"coveredbases"}}[$i];
		$scoregenes{$gene}{"unmaskedlength"}+=${$orderedexons{$chr}{"unmaskedbases"}}[$i];
	    }
	    else{
		$scoregenes{$gene}={"sumscore"=>${$orderedexons{$chr}{"score"}}[$i]*${$orderedexons{$chr}{"coveredbases"}}[$i],"exoniclength"=>${$orderedexons{$chr}{"coveredbases"}}[$i], "unmaskedlength"=>${$orderedexons{$chr}{"unmaskedbases"}}[$i]};	
	    }
	}
	
	
	foreach my $gene (keys %scoregenes){
	    my $meanscore="NA";
	    
	    if($scoregenes{$gene}{"exoniclength"}>0){
		$meanscore=($scoregenes{$gene}{"sumscore"}+0.0)/($scoregenes{$gene}{"exoniclength"}+0.0);
	    }
	    
	    my $start=min @{$exons{$gene}{"start"}};
	    my $end=max @{$exons{$gene}{"end"}};
	    my $strand=$exons{$gene}{"strand"};
	    
	    print $outputgene $gene."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$meanscore."\t".$scoregenes{$gene}{"exoniclength"}."\t".$scoregenes{$gene}{"unmaskedlength"}."\n";
	}
	
    }
}

close($outputex);
close($outputgene);
print "Done.\n";


##############################################################
