use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX;
use strict;

##############################################################

sub readCoveredRegions{

    my $pathin=$_[0];
    my $scores=$_[1];

    open(my $input, $pathin);
    
    my %header;
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $chr=$s[$header{"Chr"}];
 	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	
	if($start!=0 && $end!=0){
	    if(!(exists $scores->{$chr})){
		$scores->{$chr}={};
	    }
	    
	    if(!(exists $scores->{$chr}{$strand})){
		$scores->{$chr}{$strand}={"start"=>[], "end"=>[]};
	    }
	    
	    push(@{$scores->{$chr}{$strand}{"start"}}, $start);
	    push(@{$scores->{$chr}{$strand}{"end"}}, $end);
	    
	}

	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub combineRegions{
    my $unordered=$_[0];
    my $ordered=$_[1];

    foreach my $chr (keys %{$unordered}){
	$ordered->{$chr}={};

	foreach my $strand (keys %{$unordered->{$chr}}){
	    $ordered->{$chr}{$strand}={"start"=>[], "end"=>[]};

	    my %hashpos;
	    my $nbstart=@{$unordered->{$chr}{$strand}{"start"}};

	    for(my $i=0; $i<$nbstart; $i++){
		my $start=${$unordered->{$chr}{$strand}{"start"}}[$i];
		my $end=${$unordered->{$chr}{$strand}{"end"}}[$i];

		if(exists $hashpos{$start}){
		    if($end > $hashpos{$start}){
			$hashpos{$start}=$end;
		    }
		}
		else{
		    $hashpos{$start}=$end;
		}
	    }

	    my @uniquestart=keys %hashpos;
	    my @sortedstart=sort {$a<=>$b} @uniquestart;

	    my $currentstart=$sortedstart[0];
	    my $currentend=$hashpos{$currentstart};
	    my $nbpos=@sortedstart;

	    for(my $i=1; $i<$nbpos; $i++){
		my $thisstart=$sortedstart[$i];
		my $thisend=$hashpos{$thisstart};

		if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		    if($thisend>$currentend){
			$currentend=$thisend;
		    }
		}
		else{
		    push(@{$ordered->{$chr}{$strand}{"start"}}, $currentstart);
		    push(@{$ordered->{$chr}{$strand}{"end"}}, $currentend);
		    
		    $currentstart=$thisstart;
		    $currentend=$thisend;
		}
	    }
	   
	    ## don't forget last element
	    
	    push(@{$ordered->{$chr}{$strand}{"start"}}, $currentstart);
	    push(@{$ordered->{$chr}{$strand}{"end"}}, $currentend);
	}
    }
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script combines CSF-covered scores.\n";
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

$parameters{"pathsRegions"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathsRegions", "pathOutput");


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

print "Reading CSF-covered regions..\n";
my %scores;

my @paths=split(",", $parameters{"pathsRegions"});

foreach my $path (@paths){
    print $path."\n";
 
    readCoveredRegions($path, \%scores);
}

print "Done.\n";

##############################################################

print "Combining regions...\n";

my %orderedscores;

combineRegions(\%scores, \%orderedscores);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "Chr\tStart\tEnd\tStrand\n";

foreach my $chr (keys %orderedscores){
    foreach my $strand (keys %{$orderedscores{$chr}}){
	my $nb=@{$orderedscores{$chr}{$strand}{"start"}};
	
	for(my $i=0; $i<$nb; $i++){
	    my $start=${$orderedscores{$chr}{$strand}{"start"}}[$i];
	    my $end=${$orderedscores{$chr}{$strand}{"end"}}[$i];

	    print $output $chr."\t".$start."\t".$end."\t".$strand."\n";
	}
    }
}

close($output);

print "Done.\n";

##############################################################
