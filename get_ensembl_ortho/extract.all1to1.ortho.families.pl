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

	if(exists $splist->{$sp}){
	    $idsp->{$prefix}=$sp;
	}
	
	$line=<$input>;
    }
    
    close($input);
}

##########################################################################

sub readOrthoPairs{
    my $pathin=$_[0];
    my $idsp=$_[1];
    my $connections=$_[2];

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
	    $nbkept++;
	    my @ids=@{$ortho{$idhom}};
	    my $id1=$ids[0];
	    my $id2=$ids[1];

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

	delete $ortho{$idhom};
    }

    print "We kept ".$nbkept." orthology relationships between the interesting species.\n";
    
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

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This reconstructs families of 1-1 orthologues, using information from Ensembl.\n";
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
$parameters{"speciesList"}="NA";
$parameters{"prefixGeneInfo"}="NA";
$parameters{"suffixGeneInfo"}="NA";
$parameters{"pathOutputFamilies"}="NA";

my %defaultvalues;
my @defaultpars=("pathHomology","pathEnsemblIDs","speciesList","prefixGeneInfo","suffixGeneInfo","pathOutputFamilies");

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

my @all=split(",",$parameters{"speciesList"});
my $nball=@all;

my %allsp;

foreach my $sp (@all){
    $allsp{$sp}=1;
}
print "There are ".$nball." species in total.\n";

print "Reading gene biotypes for each species...\n";
my %biotypes;

foreach my $sp (keys %allsp){
    $biotypes{$sp}={};

    my $path=$parameters{"prefixGeneInfo"}.$sp."/".$parameters{"suffixGeneInfo"};

    if(-e $path){
	readGeneBiotype($path,$biotypes{$sp});
    }
    else{
	print "cannot find file: ".$path."\n";
    }
}

print "Done.\n";

#####################################################################

print "Reading species IDs...\n";
my %idsp;
readSpeciesIDs($parameters{"pathEnsemblIDs"},\%allsp,\%idsp);
print "Done.\n\n";

print "Reading homology relationships...\n";
my %connections;
readOrthoPairs($parameters{"pathHomology"},\%idsp,\%connections);
print "Done.\n\n";

print "Constructing clusters...\n";
my %clusters;
my %clusterids;
extractClusters(\%connections,\%clusters,\%clusterids);
print "Done.\n\n";

#####################################################################

print "Writing output for 1-1 families...\n";

my $nbnot1to1=0;
my $nbnotallpresent=0;
my $nbnotallpc=0;

my %nbspecies;

open(my $output,">".$parameters{"pathOutputFamilies"});

my @species=keys %allsp;

my $line="";

foreach my $sp (@species){
    $line.=$sp."\t";
}

chop $line;

print $output $line."\n";

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

	my $allpresent=1;
	


	foreach my $sp (keys %allsp){
	    if(!(exists $genes{$sp})){
		$allpresent=0;
	    }
	}
	
	if($allpresent==1){
	    
	    my $allpc=1;
	    
	    foreach my $sp (keys %allsp){
		my $id=$genes{$sp};
		
		if(exists $biotypes{$sp}{$id}){
		    my $bio=$biotypes{$sp}{$id};
		    
		    if($bio ne "protein_coding"){
			$allpc=0;
			last;
		    }
		}
		else{
		    print "Cannot find ".$id." in gene info for ".$sp."\n";
		}
	    }
	    
	    if($allpc==1){
		my $line="";
		
		foreach my $sp (@species){
		    $line.=$genes{$sp}."\t";
		}
		
		chop $line;
		print $output $line."\n";
	    }
	    else{
		$nbnotallpc++;
	    }
	    
	}
	else{
	    my $nbpres=keys %genes;
	    
	    if(exists $nbspecies{$nbpres}){
		$nbspecies{$nbpres}++;
	    }
	    else{
		$nbspecies{$nbpres}=1;
	    }
	    
	    $nbnotallpresent++;
	}
    }
    else{
	$nbnot1to1++;
    }
}

print "We rejected ".$nbnot1to1." families that were not all 1-1.\n";
print "We rejected ".$nbnotallpresent." families that did not have all species.\n";
print "We rejected ".$nbnotallpc." families that were not annotated as protein-coding genes.\n";

my $nballsp=@species;

for(my $i=2; $i<$nballsp; $i++){
    print "There are ".$nbspecies{$i}." families with ".$i." species.\n";
}

close($output);
print "Done.\n\n";
