#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readSwissProt{
    my $path=$_[0];
    my $acceptedscores=$_[1];
    my $reffasta=$_[2];
    
    my %rejected;
    my $nbna=0;
   
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

	    my $pescore=findInfo("PE\=",\@s);

	    if($pescore eq "NA"){
		print $id." ".$pescore."\n";
		$nbna++;
	    }

	    if(!exists $acceptedscores->{$pescore}){
		$rejected{$id}=1;
	    }

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

    print "Found ".$nbna." proteins with no PE score.\n";

    my $nbrejected=keys %rejected;

    foreach my $id (keys %rejected){
	delete $reffasta->{$id};
    }
    
    print "Removed ".$nbrejected." proteins with bad PE scores.\n";
}


###############################################################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];
    
    my $res="NA";
    
    my @grepres=grep(/${pattern}/,@{$array});
    
    if(@grepres==1){
	my @t=split("\=",$grepres[0]);
	$res=$t[1];
    }
    
    return $res;
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

###############################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script filters SwissProt proteins by PE score. \n";
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

$parameters{"pathSwissProt"}="NA";
$parameters{"acceptedScores"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathSwissProt", "acceptedScores", "pathOutput");

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

print "Extracting acceptable scores...\n";

my %scores;
my @s=split(",", $parameters{"acceptedScores"});

for(my $i=0; $i<@s; $i++){
    $scores{$s[$i]}=1;
}

my $nbscores=keys %scores;

print "There are ".$nbscores." accepted scores: ".join("; ",keys %scores)."\n";

print "Done.\n";

#####################################################################################

print "Reading and filtering SwissProt sequences...\n";

my %swissprot;

readSwissProt($parameters{"pathSwissProt"}, \%scores, \%swissprot);

my $nbprot=keys %swissprot;

print "Kept ".$nbprot." sequences.\n";
print "Done.\n";

#####################################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

foreach my $id (keys %swissprot){
    writeSequence($swissprot{$id}, $id, $output);
}

close($output);

print "Done.\n";

#####################################################################################
