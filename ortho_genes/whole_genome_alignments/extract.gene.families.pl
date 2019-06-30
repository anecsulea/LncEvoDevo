use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###########################################################################
###########################################################################

sub readConnectedGenes{
    my $pathin=$_[0];
    my $connectedgenes=$_[1];
    my $allspecies=$_[2];

    open(my $input, $pathin);
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my @species;

    for(my $i=0; $i<@s; $i++){
	my @t=split("\\.",$s[$i]);

	if($t[0] ne "ID"){
	    print "Weird header: ".$line."\n";
	}

	my $sp=$t[1];

	push(@species, $sp);

	if(!exists $allspecies->{$sp}){
	    $allspecies->{$sp}=1;
	}
    }

    $line=<$input>;

    my $nbsp=@species;

    print "Found ".$nbsp." species in file: ".join(",", @species)."\n";

    if($nbsp<2){
	print "Weird number of species!\n";
	exit(1);
    }
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my @ids;

	for(my $i=0; $i<$nbsp; $i++){
	    my $fullid=$species[$i].":".$s[$i];
	    push(@ids, $fullid);
	}

	for(my $i=0; $i<$nbsp; $i++){
	    my $id1=$ids[$i];

	    if(!exists $connectedgenes->{$id1}){
		$connectedgenes->{$id1}={};
	    }
	    
	    for(my $j=0; $j<$nbsp; $j++){
		my $id2=$ids[$j];
		$connectedgenes->{$id1}{$id2}=1;
	    }
	}

	$line=<$input>;
    }

    close($input);
}

###########################################################################

sub addToCluster{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $indexclusters=$_[2];
    my $key=$_[3];

    if(exists $refconnected->{$key}){
                    
        my $indexcluster="NA";
	
        if(exists $indexclusters->{$key}){
            $indexcluster=$indexclusters->{$key}
        }       
              
        if($indexcluster eq "NA"){
            my $nbclusters=keys %{$refclusters};
            $indexcluster=$nbclusters+1;
            $refclusters->{$indexcluster}=[$key];
            $indexclusters->{$key}=$indexcluster;
        }
                
        foreach my $connection (keys %{$refconnected->{$key}}){
            if(!(exists $indexclusters->{$connection})){
                push(@{$refclusters->{$indexcluster}},$connection);
                $indexclusters->{$connection}=$indexcluster;
                addToCluster($refconnected,$refclusters,$indexclusters,$connection);
            }
        }

        delete $refconnected->{$key};
    }
}

###########################################################################

sub extractClusters{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    
    my %indexclusters;
    
    my $nbconnected=keys %{$refconnected};
    
    my $round=0;
    
    my %alreadyprinted;
    
    while($nbconnected>0){
	
        foreach my $key (keys %{$refconnected}){
            addToCluster($refconnected,$refclusters,\%indexclusters,$key);
            $nbconnected=keys %{$refconnected};
        }
        
        $round++;
        
        $nbconnected=keys %{$refconnected};
    }
}

###########################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts gene families starting from reciprocal best hits.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

###########################################################################
###########################################################################

my %parameters;

$parameters{"pathsReciprocalBestHits"}="NA";
$parameters{"pathOutput1to1Families"}="NA";
$parameters{"pathOutputMultipleFamilies"}="NA";

my @defaultpars=("pathsReciprocalBestHits", "pathOutput1to1Families", "pathOutputMultipleFamilies");

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

print "Reading connected genes...\n";

my @paths=split(",", $parameters{"pathsReciprocalBestHits"});

my %connectedgenes;
my %allspecies;

foreach my $path (@paths){
    if($path ne ""){
	print "Reading from ".$path."\n";
	
	readConnectedGenes($path, \%connectedgenes, \%allspecies);
	my $nbc=keys %connectedgenes;
	
	print $nbc." connected genes so far.\n";
    }
}

my @species=keys %allspecies;
my $nbsp=@species;

print "Found ".$nbsp." species: ".join(",", @species)."\n";

print "Done.\n";

##############################################################

print "Extracting clusters...\n";

my %geneclusters;
extractClusters(\%connectedgenes, \%geneclusters);

print "Done.\n";

##############################################################

print "Analyzing clusters and writing output...\n";

open(my $output1to1, ">".$parameters{"pathOutput1to1Families"});
open(my $outputmulti, ">".$parameters{"pathOutputMultipleFamilies"});

my $header="";

foreach my $sp (@species){
    $header.="ID.".$sp."\t";
}

chop $header;

print $output1to1 $header."\n";
print $outputmulti $header."\n";

foreach my $indexclust (keys %geneclusters){
    my %genessp;

    my $multi=0;

    foreach my $gene (@{$geneclusters{$indexclust}}){
	my @s=split(":", $gene);
	my $sp=$s[0];
	my $g=$s[1];

	if(exists $genessp{$sp}){
	    push(@{$genessp{$sp}}, $g);
	    $multi=1;
	} else{
	    $genessp{$sp}=[$g];
	}
    }

    
    my $line="";
    
    foreach my $sp (@species){
	if(exists $genessp{$sp}){
	    $line.=join(";",@{$genessp{$sp}})."\t";
	} else{
	    $line.="NA\t";
	}
    }
    
    chop $line;

    if($multi==1){
	print $outputmulti $line."\n";
    } else{
	print $output1to1 $line."\n";
    }
}

close($output1to1);
close($outputmulti);

print "Done.\n";

##############################################################
