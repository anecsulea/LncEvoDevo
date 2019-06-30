use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub extractClusters{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $refclustid=$_[2];

    my $nbconnected=keys %{$refconnected};

    my $round=0;

    while($nbconnected>0){
	
	foreach my $key (keys %{$refconnected}){
	    addToCluster($refconnected,$refclusters,$refclustid,$key);
	}
	
	$round++;
	
	$nbconnected=keys %{$refconnected};
    }
}

##############################################################

sub addToCluster{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $refclustid=$_[2];
    my $key=$_[3];

    if(exists $refconnected->{$key}){
	
	## find the cluster that contains this key
	
	my $indexcluster="NA";
	
	if(exists $refclustid->{$key}){
	    $indexcluster=$refclustid->{$key};
	}
	
	## if there isn't any
	
	if($indexcluster eq "NA"){
	    my $nbclusters=keys %{$refclusters};
	    $indexcluster=$nbclusters+1;
	    $refclusters->{$indexcluster}={$key=>1};
	    $refclustid->{$key}=$indexcluster;
	}
		
	foreach my $connection (keys %{$refconnected->{$key}}){

	    ## check if this island is already in the cluster
	    
	    if(!(exists $refclusters->{$indexcluster}{$connection})){
		$refclusters->{$indexcluster}{$connection}=1;
		$refclustid->{$connection}=$indexcluster;
		addToCluster($refconnected,$refclusters,$refclustid,$connection);
	    }
	}

	## after we've checked all of its connections, remove it from the connected islands

	delete $refconnected->{$key};
    }
    
}

##########################################################################

sub readSpeciesIDs{
    my $pathin=$_[0];
    my $splist=$_[1];
    my $idsp=$_[2];
    
    open(my $input, $pathin);

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split(" ",$line);
	my $sp=$s[0];
	my $prefix=$s[1];
	
	$idsp->{$prefix}=$sp;
	
	$splist->{$sp}=1;
	
	$line=<$input>;
    }
    
    close($input);
}

##########################################################################

sub readOrthoPairs{
    my $pathin=$_[0];
    my $idsp=$_[1];
    my $existingfam=$_[2];
    my $connections=$_[3];

    my %ortho;

    open(my $input, $pathin);

    my $line=<$input>; ## header
    $line=<$input>;
    
    my $nbdone=0;

    while($line){
	chomp $line;
	my @s=split("\t",$line);
	my $gene=$s[0];

	if($gene =~ m{\d+}){
	    my $posdigit=$-[0];
	    my $prefix=substr $gene,0,$posdigit;
	   	    
	    if(exists $idsp->{$prefix}){
		my $idhom=$s[1];
		my $sp=$idsp->{$prefix};
		my $id=$sp.":".$gene;
		
		if(exists $ortho{$idhom}){
		    push(@{$ortho{$idhom}},$id);
		}
		else{
		    $ortho{$idhom}=[$id];
		}
	    }
	}
	else{
	    print $gene. "has no digits";
	}
	
	$nbdone++;

	if($nbdone%1000000==0){
	    print $nbdone. " lines read\n";
	}

	$line=<$input>;
    }

    my $nbhom=keys %ortho;
    print "There were ".$nbhom." orthology relationships initially.\n";
    my $nbkept=0;
    
    foreach my $idhom (keys %ortho){
	my $nb=@{$ortho{$idhom}};

	if($nb==2){

	    my @ids=@{$ortho{$idhom}};
	    my $id1=$ids[0];
	    my $id2=$ids[1];

	    my @s1=split(":",$id1);
	    my @s2=split(":",$id2);
	    
	    if(exists $existingfam->{$s1[1]} || exists $existingfam->{$s2[1]}){
		$nbkept++;
		
		if(exists $connections->{$id1}){
		    $connections->{$id1}{$id2}=1;
		}
		else{
		    $connections->{$id1}={$id2=>1};
		}
		
		if(exists $connections->{$id2}){
		    $connections->{$id2}{$id1}=1;
		}
		else{
		    $connections->{$id2}={$id1=>1};
		}
	    }
	}

	delete $ortho{$idhom};
    }

    print "We kept ".$nbkept." orthology relationships between the interesting species, that connect with the existing families.\n";
    
    close($input);
}

##############################################################

sub readGeneBiotype{

    my $pathBiotype=$_[0];
    my $bio=$_[1];
    
    open(my $input,$pathBiotype);
    
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

	$bio->{$id}=$thisbio;

	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub readCoreFamilies{
    my $pathin=$_[0];
    my $families=$_[1];
    my $idfamilies=$_[2];

    open(my $input, $pathin);
    
    my $line=<$input>; ## header
    $line=<$input>; ## header
    my $index=0;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	$families->{$index}=[];
	push(@{$families->{$index}}, @s);
	
	foreach my $id (@s){
	    $idfamilies->{$id}=1;
	}

	$line=<$input>;
	$index++;
    }

    close($input);
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This reconstructs adds outgroup species to families of 1-1 orthologues.\n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################
##########################################################################

## define parameters

my %parameters;
$parameters{"pathHomology"}="NA";
$parameters{"pathEnsemblIDs"}="NA";
$parameters{"pathCoreFamilies"}="NA";
$parameters{"prefixGeneInfo"}="NA";
$parameters{"suffixGeneInfo"}="NA";
$parameters{"firstIndex"}=0;
$parameters{"pathOutputFamilies"}="NA";

my %defaultvalues;
my @defaultpars=("pathHomology","pathEnsemblIDs","pathCoreFamilies","prefixGeneInfo","suffixGeneInfo","firstIndex","pathOutputFamilies");

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
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

print "Reading species IDs...\n";
my %idsp;
my %allsp;
readSpeciesIDs($parameters{"pathEnsemblIDs"},\%allsp,\%idsp);
my $nbsp=keys %allsp;
print "There are ".$nbsp." species in total.\n";
print "Done.\n\n";


print "Reading gene biotypes for each species...\n";
my %biotypes;

foreach my $sp (keys %allsp){
    $biotypes{$sp}={};

    my $path=$parameters{"prefixGeneInfo"}.$sp."/".$parameters{"suffixGeneInfo"};

    if(-e $path){
	readGeneBiotype($path,$biotypes{$sp});
    }
    else{
	print "cannot find gene info file for ".$sp."\n";
    }
}
print "Done.\n";

#####################################################################

print "Reading core families...\n";
my %existingfamilies;
my %idfamilies;
readCoreFamilies($parameters{"pathCoreFamilies"}, \%existingfamilies,\%idfamilies);
my $nbfam=keys %existingfamilies;
print "There are ".$nbfam." core families.\n";
print "Done.\n";

print "Reading homology relationships...\n";
my %connections;
readOrthoPairs($parameters{"pathHomology"},\%idsp,\%idfamilies,\%connections);
print "Done.\n\n";

print "Constructing clusters...\n";
my %clusters;
my %clusterids;
extractClusters(\%connections,\%clusters,\%clusterids);
print "Done.\n\n";

#####################################################################

print "Writing output for 1-1 families...\n";

my $firstIndex=$parameters{"firstIndex"}+0;
print "First family will be named Family".$firstIndex."\n";

my $nbnot1to1=0;
my $nbnotallpc=0;

my %nbspecies;

open(my $output,">".$parameters{"pathOutputFamilies"});

my @species=keys %allsp;

my $line="FamilyID";

foreach my $sp (@species){
    $line.="\t".$sp;
}

print $output $line."\n";

my $index=$firstIndex;

foreach my $clustid (keys %clusters){
    my %genes;
    my $all1to1=1;
      
    foreach my $gene (keys  %{$clusters{$clustid}}){
	my @s=split(":",$gene);
	my $sp=$s[0];
	my $id=$s[1];

	if(exists $genes{$sp}){
	    $all1to1=0;
	    last;
	}
	else{
	    $genes{$sp}=$id;
	}
    }

    if($all1to1==1){

	my $allpc=1;
	
	foreach my $sp (keys %genes){
	    if(exists $biotypes{$sp}){
	       my $id=$genes{$sp};
	       
	       if(exists $biotypes{$sp}{$id}){
		   my $bio=$biotypes{$sp}{$id};
		
		   if($bio ne "protein_coding"){
		       $allpc=0;
		    last;
		   }
	       }
	       # else{
	       # 	   print "Cannot find ".$id." in gene info for ".$sp."\n";
	       # }
	    }
	}
	
	if($allpc==1){
	    my $famid="Family".$index;

	    my $line=$famid;
	    
	    foreach my $sp (@species){
		if(exists $genes{$sp}){
		    $line.="\t".$genes{$sp};
		}
		else{
		    $line.="\tNA";
		}
	    }
	    
	    print $output $line."\n";

	    $index++;
	}
	else{
	    $nbnotallpc++;
	}
    }
    else{
	$nbnot1to1++;
    }
}

print "We rejected ".$nbnot1to1." families that were not all 1-1.\n";
print "We rejected ".$nbnotallpc." families that were not annotated as protein-coding genes.\n";


close($output);
print "Done.\n\n";
