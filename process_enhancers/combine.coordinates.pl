use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readEnhancers{
    my $pathin=$_[0];
    my $halfsize=$_[1];
    my $coords=$_[2];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];
    my $input;
    
    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[0];
	my $midpoint=$s[1];
	my $start=$midpoint-$halfsize;
	my $end=$midpoint+$halfsize;
	
	if(exists $coords->{$chr}){
	    if(exists $coords->{$chr}{$start}){
		if($end>$coords->{$chr}{$start}){
		    $coords->{$chr}{$start}=$end;
		}
	    } else{
		$coords->{$chr}{$start}=$end;
	    }
	}
	else{
	    $coords->{$chr}={$start=>$end};
	}
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub combineEnhancers{
    my $coords=$_[0];
    my $combined=$_[1];

    foreach my $chr (keys %{$coords}){
	$combined->{$chr}={"start"=>[], "end"=>[]};

	my @uniquestart=keys %{$coords->{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	my $nbstart=@sortedstart;
	my $currentstart=$sortedstart[0];
	my $currentend=$coords->{$chr}{$currentstart};

	for(my $i=1; $i<$nbstart; $i++){
	    my $thisstart=$sortedstart[$i];
	    my $thisend=$coords->{$chr}{$thisstart};

	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    } else{
		push(@{$combined->{$chr}{"start"}}, $currentstart);
		push(@{$combined->{$chr}{"end"}}, $currentend);

		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}

	push(@{$combined->{$chr}{"start"}}, $currentstart);
	push(@{$combined->{$chr}{"end"}}, $currentend);
    }
}

##############################################################

sub printHelp{
    
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script combines enhancer coordinates.\n";
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
$parameters{"pathsEnhancers"}="NA";
$parameters{"halfSize"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathsEnhancers","halfSize","pathOutput");
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

print "Reading enhancer coordinates...\n";

my $halfsize=$parameters{"halfSize"}+0;
print "half size: ".$halfsize."\n";

my %enhancers;
my @paths=split(",", $parameters{"pathsEnhancers"});

foreach my $path (@paths){
    if(-e $path){
	print "reading from ".$path."\n";
	readEnhancers($path, $halfsize, \%enhancers);
    }
}

print "Done.\n";

print "Combining enhancers...\n";
my %combined;
combineEnhancers(\%enhancers, \%combined);
print "Done.\n";

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});
foreach my $chr (keys %combined){
    my $nbe=@{$combined{$chr}{"start"}};
    
    for(my $i=0; $i<$nbe; $i++){
	print $output $chr."\t".${$combined{$chr}{"start"}}[$i]."\t".${$combined{$chr}{"end"}}[$i]."\n";
    }
}

close($output);

print "Done.\n";

##############################################################
##############################################################
