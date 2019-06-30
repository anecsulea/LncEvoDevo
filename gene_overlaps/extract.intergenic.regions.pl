use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readBlocks{
    my $pathin=$_[0];
    my $blocks=$_[1];

    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $gene=$s[0];
	my $chr=$s[2];
	my $start=$s[3];
	my $end=$s[4];
	my $strand=$s[5];

	if(!(exists $blocks->{$gene})){
	    $blocks->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>[], "end"=>[]};
	}

	push(@{$blocks->{$gene}{"start"}}, $start);
	push(@{$blocks->{$gene}{"end"}}, $end);
	
	$line=<$input>;
    }

    close($input);
}

###################################################################

sub extractGeneCoords{
    my $blocks=$_[0];
    my $genecoords=$_[1];

    foreach my $gene (keys %{$blocks}){
	my $chr=$blocks->{$gene}{"chr"};
	my $strand=$blocks->{$gene}{"strand"};

	my $start=min @{$blocks->{$gene}{"start"}};
	my $end=max @{$blocks->{$gene}{"end"}};

	$genecoords->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>$start, "end"=>$end};
    }
}

###########################################################################

sub makeGeneBlocks{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $gene (keys %{$coords}){
	my $chr=$coords->{$gene}{"chr"};
	my $strand=$coords->{$gene}{"strand"};
	my $start=$coords->{$gene}{"start"};
	my $end=$coords->{$gene}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		if($end>$hashpos{$chr}{$start}){
		    $hashpos{$chr}{$start}=$end;
		}
	    }
	    else{
		$hashpos{$chr}{$start}=$end;
	    }
	}
	else{
	    $hashpos{$chr}={$start=>$end};
	}
    }
   

    foreach my $chr (keys %hashpos){
	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	$ordered->{$chr}={"start"=>[],"end"=>[]};

	my $nbpos=@sortedstart;
	my $currentstart=$sortedstart[0];
	my $currentend=$hashpos{$chr}{$currentstart};

	for(my $i=1; $i<$nbpos; $i++){
	    my $thisstart=$sortedstart[$i];
	    my $thisend=$hashpos{$chr}{$thisstart};

	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    } else{
		push(@{$ordered->{$chr}{"start"}}, $currentstart);
		push(@{$ordered->{$chr}{"end"}}, $currentend);

		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}

	push(@{$ordered->{$chr}{"start"}}, $currentstart);
	push(@{$ordered->{$chr}{"end"}}, $currentend);

    }
}

###################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts intergenic region coordinates. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

###################################################################

my %parameters;

$parameters{"pathExonBlocks"}="NA";
$parameters{"minDistance"}="NA";
$parameters{"minSize"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonBlocks","minDistance", "minSize", "pathOutput");

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

print "Reading exon block coordinates for all genes...\n";

my %allblocks;

readBlocks($parameters{"pathExonBlocks"}, \%allblocks);

my $nbg=keys %allblocks;

print "There are ".$nbg." genes in total.\n";

my %genecoords;

extractGeneCoords(\%allblocks, \%genecoords);  

print "Done.\n";

#####################################################################

print "Constructing gene blocks...\n";
my %geneblocks;
makeGeneBlocks(\%genecoords, \%geneblocks);
print "Done.\n";

#####################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

my $mindist=$parameters{"minDistance"}+0;
my $minsize=$parameters{"minSize"}+0;
print "minimum distance: ".$mindist."\n";
print "minimum size: ".$minsize."\n";

foreach my $chr (keys %geneblocks){
    my $nbg=@{$geneblocks{$chr}{"start"}};

    if($nbg>=2){
	for(my $i=0; $i<($nbg-1); $i++){
	    my $end1=${$geneblocks{$chr}{"end"}}[$i];
	    my $start2=${$geneblocks{$chr}{"start"}}[$i+1];

	    if($end1>=$start2){
		print "Weird! these are not gene blocks: ".$end1." ".$start2."\n";
		exit(1);
	    }

	    my $size=$start2-$end1-1;

	    if($size>=(2*$mindist+$minsize)){
		my $startreg=$end1+$mindist;
		my $endreg=$start2-$mindist;

		print $output $chr."\t".$startreg."\t".$endreg."\n";
	    }
	    
	}
    }
}

close($output);

print "Done.\n";

#####################################################################
