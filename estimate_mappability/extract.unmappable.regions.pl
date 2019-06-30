use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readMappedRegions{
    my $path=$_[0];
    my $map=$_[1];

    open(my $input, $path);

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);
	my $chr=$s[0];
	my $start=$s[1]+0;
	my $end=$s[2]+0;
	
	if(exists $map->{$chr}){
	    push(@{$map->{$chr}{"start"}},$start);
	    push(@{$map->{$chr}{"end"}},$end);
	}
	else{
	    $map->{$chr}={"start"=>[$start],"end"=>[$end]};
	}
	
	$line=<$input>;
    }
    close($input);
}

##############################################################

sub orderRegions{
    my $unordered=$_[0];
    my $ordered=$_[1];

    foreach my $chr (keys %{$unordered}){
	$ordered->{$chr}={"start"=>[],"end"=>[]};

	my %hashpos;

	my $nbpos=@{$unordered->{$chr}{"start"}};

	for(my $i=0; $i<$nbpos; $i++){
	    my $start=${$unordered->{$chr}{"start"}}[$i];
	    my $end=${$unordered->{$chr}{"end"}}[$i];
	    
	    $hashpos{$start}=$end;
	}

	my @uniquestart=keys %hashpos;
	my @sortedstart= sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    push(@{$ordered->{$chr}{"start"}},$start);
	    push(@{$ordered->{$chr}{"end"}},$hashpos{$start});
	}
    }
}

##############################################################

sub extractUnmappedRegions{
    my $mapped=$_[0];
    my $chrsizes=$_[1];
    my $unmapped=$_[2];
    
    foreach my $chr (keys %{$mapped}){

	my $chrsize=$chrsizes->{$chr};

	$unmapped->{"chr"}={"start"=>[],"end"=>[]};
    
	my $nbmap=@{$mapped->{$chr}{"start"}};

	my $firststart=${$mapped->{$chr}{"start"}}[0];
	
	if($firststart>1){
	    push(@{$unmapped->{$chr}{"start"}},1);
	    push(@{$unmapped->{$chr}{"end"}},($firststart-1));
	}
	
	for(my $i=0; $i<($nbmap-1); $i++){
	    my $thisend=${$mapped->{$chr}{"end"}}[$i];
	    my $nextstart=${$mapped->{$chr}{"start"}}[$i+1];
	    
	    if(($thisend+1)<$nextstart){
		push(@{$unmapped->{$chr}{"start"}},$thisend+1);
		push(@{$unmapped->{$chr}{"end"}},$nextstart-1);
	    }
	    else{
		print "Weird!!! we should have collapsed this entry ".$thisend." ".$nextstart."\n";
	    }
	}
	
	my $lastend=${$mapped->{$chr}{"end"}}[-1];
	
	if($lastend<$chrsize){
	    push(@{$unmapped->{$chr}{"start"}},$lastend+1);
	    push(@{$unmapped->{$chr}{"end"}},$chrsize);
	}
    }
}

##############################################################

sub readChromosomeSizes{
    my $pathin=$_[0];
    my $chrsizes=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);
	my $chr=$s[0];
	my $size=$s[1]+0;
	$chrsizes->{$chr}=$size;
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub printHelp{
    
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts unmappable regions of the genome.\n";
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
$parameters{"pathsMappedRegions"}="NA";
$parameters{"pathChromosomeSizes"}="NA";
$parameters{"pathOutputTxt"}="NA";
$parameters{"pathOutputBedGraph"}="NA";

my @defaultpars=("pathsMappedRegions","pathChromosomeSizes","pathOutputTxt","pathOutputBedGraph");
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

print "Reading mapped regions...\n";
my %unorderedmap;
my @paths=split(",",$parameters{"pathsMappedRegions"});
foreach my $path (@paths){
    if(-e $path){
	readMappedRegions($path,\%unorderedmap);
    }
    else{
	print "Cannot find path: ".$path."\n";
    }
}
my $nbchr=keys %unorderedmap;
print "There are mapped regions on ".$nbchr." chromosomes.\n";
print "Done.\n";

##############################################################

print "Ordering mapped regions...\n";
my %mapped;
orderRegions(\%unorderedmap, \%mapped);
print "Done.\n";

##############################################################

print "Reading chromosome sizes...\n";
my %chrsizes;
readChromosomeSizes($parameters{"pathChromosomeSizes"},\%chrsizes);
my $nbchr=keys %chrsizes;
print "Done.\n";

##############################################################

print "Extracting unmapped regions...\n";
my %unmapped;
extractUnmappedRegions(\%mapped, \%chrsizes, \%unmapped);
print "Done.\n";

print "Writing output...\n";
open(my $output, ">".$parameters{"pathOutputTxt"});
foreach my $chr (keys %unmapped){
    my $nbreg=@{$unmapped{$chr}{"start"}};
    
    for(my $i=0; $i<$nbreg; $i++){
	my $start=${$unmapped{$chr}{"start"}}[$i];
	my $end=${$unmapped{$chr}{"end"}}[$i];
	print $output $chr."\t".$start."\t".$end."\n";
    }
}
close($output);
print "Done.\n";


print "Writing output...\n";
open(my $output, ">".$parameters{"pathOutputBedGraph"});
foreach my $chr (keys %unmapped){
    my $nbreg=@{$unmapped{$chr}{"start"}};
    
    for(my $i=0; $i<$nbreg; $i++){
	my $start=${$unmapped{$chr}{"start"}}[$i];
	my $end=${$unmapped{$chr}{"end"}}[$i];
	print $output $chr."\t".($start-1)."\t".$end."\n";
    }
}
close($output);
print "Done.\n";

##############################################################
