use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readIntergenicRegions{
    my $pathin=$_[0];
    my $regions=$_[1];

    $regions->{"chr"}=[];
    $regions->{"start"}=[];
    $regions->{"end"}=[];
    
    open(my $input, $pathin);

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[0];
	my $start=$s[1]+0;
	my $end=$s[2]+0;

	push(@{$regions->{"chr"}}, $chr);
	push(@{$regions->{"start"}}, $start);
	push(@{$regions->{"end"}}, $end);

	$line=<$input>;
    }
    
    close($input);
}

############################################################################

sub sampleRegions{
    my $regions=$_[0];
    my $size=$_[1];
    my $flank=$_[2];
    my $pathout=$_[3];

    open(my $output, ">".$pathout);

    my $nbreg=@{$regions->{"chr"}};

    my $minsize=$size+2*$flank;
    
    for(my $i=0; $i<$nbreg; $i++){
	my $chr=${$regions->{"chr"}}[$i];
	my $rstart=${$regions->{"start"}}[$i];
	my $rend=${$regions->{"end"}}[$i];

	my $totsize=$rend-$rstart+1;

	if($totsize>=$minsize){
	    my $startpos=$rstart+int(rand($totsize-$size-$flank));
	    my $endpos=$startpos+$size-1;

	    print $output $chr."\t".$startpos."\t".$endpos."\n";
	}
    }
        
    close($output);
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

$parameters{"pathIntergenicRegions"}="NA";
$parameters{"regionSize"}="NA";
$parameters{"minFlankingDistance"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathIntergenicRegions", "regionSize", "minFlankingDistance", "pathOutput");

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

print "Reading intergenic regions...\n";

my %igregions;
readIntergenicRegions($parameters{"pathIntergenicRegions"}, \%igregions);
my $nb=@{$igregions{"chr"}};

print "Found ".$nb." intergenic regions.\n";

print "Done.\n";

#####################################################################

my $size=$parameters{"regionSize"}+0;
my $flank=$parameters{"minFlankingDistance"}+0;

print "Resampling regions, size ".$size.", minimum flanking distance ".$flank."\n";

sampleRegions(\%igregions, $size, $flank, $parameters{"pathOutput"});

print "Done.\n";

#####################################################################
