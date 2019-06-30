use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $genecoords=$_[1];
    
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
	
	if(exists $genecoords->{$geneid}){
	    if($genecoords->{$geneid}{"start"}>$start){
		$genecoords->{$geneid}{"start"}=$start;
	    }
	    
	    if($genecoords->{$geneid}{"end"}<$end){
		$genecoords->{$geneid}{"end"}=$end;
	    }
	}
	else{
	    $genecoords->{$geneid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	}
        
	$line=<$input>;
    }

    close($input);
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

##############################################################

sub makeBlocks{
    my $intervals=$_[0];
    my $blocks=$_[1];

    foreach my $chr (keys %{$intervals}){
	my %hashpos;

	$blocks->{$chr}={"start"=>[], "end"=>[]};
	
	my $nb=@{$intervals->{$chr}{"start"}};
    
	for(my $i=0; $i<$nb; $i++){
	    my $start=${$intervals->{$chr}{"start"}}[$i];
	    my $end=${$intervals->{$chr}{"end"}}[$i];
	    
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
		push(@{$blocks->{$chr}{"start"}}, $currentstart);
		push(@{$blocks->{$chr}{"end"}}, $currentend);
		
		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}
	
	push(@{$blocks->{$chr}{"start"}}, $currentstart);
	push(@{$blocks->{$chr}{"end"}}, $currentend);
    }
}

##############################################################

sub extractIntergenicRegions{
    my $geneblocks=$_[0];
    my $igblocks=$_[1];

    foreach my $chr (keys %{$geneblocks}){
	my $nbg=@{$geneblocks->{$chr}{"start"}};
	
	if($nbg>=2){
	    $igblocks->{$chr}={"start"=>[], "end"=>[]};
	    
	    for(my $i=0; $i<($nbg-1); $i++){
		my $end=${$geneblocks->{$chr}{"end"}}[$i];
		my $start=${$geneblocks->{$chr}{"start"}}[$i+1];

		if(($end+1)<=($start-1)){
		    push(@{$igblocks->{$chr}{"start"}}, $end+1);
		    push(@{$igblocks->{$chr}{"end"}}, $start-1);
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

$parameters{"pathExonBlocks"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonBlocks", "pathOutput");

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

my %genecoords;

readExonBlocks($parameters{"pathExonBlocks"}, \%genecoords);

my $nbg=keys %genecoords;

print "Found ".$nbg." genes.\n";

print "Done.\n";

#####################################################################

print "Ordering gene coordinates...\n";

my %orderedgenes;
orderCoords(\%genecoords, \%orderedgenes);

print "Done.\n";

#####################################################################

print "Making gene blocks...\n";

my %geneblocks;
makeBlocks(\%orderedgenes, \%geneblocks);

print "Done.\n";

##############################################################

print "Extracting intergenic region coordinates...\n";

my %igregions;

extractIntergenicRegions(\%geneblocks, \%igregions);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

foreach my $chr (keys %igregions){
    my $nb=@{$igregions{$chr}{"start"}};

    for(my $i=0; $i<$nb; $i++){
	my $start=${$igregions{$chr}{"start"}}[$i];
	my $end=${$igregions{$chr}{"end"}}[$i];

	print $output $chr."\t".$start."\t".$end."\n";
    }
}

close($output);

print "Done.\n";

##############################################################
