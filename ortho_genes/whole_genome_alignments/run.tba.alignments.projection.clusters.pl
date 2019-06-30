#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
   
    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];

	    # print "saw chromosome ".$id."\n";
	    
	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

#########################################################################################

sub writeSequence{
    my $sequence=$_[0];
    my $name=$_[1];
    my $output=$_[2];

    my $n=length $sequence;

    print $output ">".$name."\n";

    my $i=0;

    while($i<($n-60)){

        my $subseq=substr $sequence,$i,60;

        print $output $subseq ."\n";

        $i+=60;
    }

    if($i<$n){
        my $subseq=substr $sequence,$i;
        print $output $subseq ."\n";
    }
}

#########################################################################################

sub readClusters{
    my $pathin=$_[0];
    my $sp1=$_[1];
    my $sp2=$_[2];
    my $clusters=$_[3];

    open(my $input, $pathin);
    my $line=<$input>; ## header
    chomp $line;

    my @s=split("\t", $line);

    my %header;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    my $indexcluster=0;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $g1=$s[$header{"Genes.".$sp1}];
	my $g2=$s[$header{"Genes.".$sp2}];

	if($g1 ne "NA" && $g2 ne "NA"){
	    my @genes1=split(";", $g1);
	    my @genes2=split(";", $g2);

	    $clusters->{$indexcluster}={$sp1=>[], $sp2=>[]};
	    push(@{$clusters->{$indexcluster}{$sp1}}, @genes1);
	    push(@{$clusters->{$indexcluster}{$sp2}}, @genes2);
	    
	    $indexcluster++;
	}
	$line=<$input>;
    }

    close($input);
}

#########################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script runs muscle alignments for projection clusters. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################################
##########################################################################################

my %parameters;

$parameters{"species1"}="NA";
$parameters{"species2"}="NA";
$parameters{"pathSequences1"}="NA";
$parameters{"pathSequences2"}="NA";
$parameters{"pathClusters"}="NA";
$parameters{"run"}="NA";
$parameters{"modulo"}="NA";
$parameters{"dirOutput"}="NA";

my %defaultvalues;
my @defaultpars=("species1", "species2", "pathSequences1", "pathSequences2", "pathClusters", "run", "modulo", "dirOutput");

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

#####################################################################################
#####################################################################################

print "Reading exon block cDNA sequences...\n";

my $sp1=$parameters{"species1"};
my $sp2=$parameters{"species2"};

my %sequences1;
my %sequences2;

readFasta($parameters{"pathSequences1"}, \%sequences1);
readFasta($parameters{"pathSequences2"}, \%sequences2);

my $nbg1=keys %sequences1;
my $nbg2=keys %sequences2;

print "Found ".$nbg1." sequences for ".$sp1."\n";
print "Found ".$nbg2." sequences for ".$sp2."\n";

print "Done.\n";

#####################################################################################

print "Reading clusters...\n";

my %clusters;

readClusters($parameters{"pathClusters"}, $sp1, $sp2, \%clusters);

my $nbclust=keys %clusters;

print "Found ".$nbclust." non-empty clusters.\n";

print "Done.\n";

#####################################################################################

print "Running TBA alignments...\n";

my $tree="(".$sp1." ".$sp2.")";

print "Tree for TBA: ".$tree."\n";


my $dirout=$parameters{"dirOutput"};

my $run=$parameters{"run"}+0;
my $modulo=$parameters{"modulo"}+0;

print "Run: ".$run." modulo ".$modulo."\n";

foreach my $indexclust (keys %clusters){
    my $remainder=0;
    
    if($modulo!=0){
    	$remainder=$indexclust%$modulo;
    }
    
    if($modulo==0 || $remainder==$run){
        foreach my $gene1 (@{$clusters{$indexclust}{$sp1}}){
	    if(!exists $sequences1{$gene1}){
		print "Weird! cannot find sequence for ".$gene1." in ".$sp1."\n";
		exit(1);
	    }
	    my $seq1=$sequences1{$gene1};
	    
	    foreach my $gene2 (@{$clusters{$indexclust}{$sp2}}){
		if(!exists $sequences2{$gene2}){
		    print "Weird! cannot find sequence for ".$gene2." in ".$sp2."\n";
		    exit(1);
		}

		my $seq2=$sequences2{$gene2};

		my $pathout=$dirout."/".$gene1.".".$gene2.".maf.gz";

		if(-e $pathout){
		    print "already done\n";
		} else{
		    my $thisoutdir=$dirout."/".$gene1.".".$gene2;
		    system("mkdir ".$thisoutdir);
		    
		    open(my $output, ">".$thisoutdir."/".$sp1);
		    writeSequence($seq1, $sp1, $output);
		    close($output);

		    open(my $output, ">".$thisoutdir."/".$sp2);
		    writeSequence($seq2, $sp2, $output);
		    close($output);

		    chdir $thisoutdir; 
		    system("all_bz \"".$tree."\"");
		    system("tba \"".$tree."\" *.*.maf tba.maf");
		    system("gzip tba.maf");
		    system("mv tba.maf.gz ".$pathout);
		    system("rm -r ".$thisoutdir);

		}
	    }
	}
    }

    if($indexclust%1000==0){
	print "Done ".$indexclust." clusters.\n";
    }
}

print "Done.\n";

#####################################################################################
