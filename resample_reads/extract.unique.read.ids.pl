use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub findNbHits{
    my $s=$_[0];

    my $nh="NA";

    my @grepres=grep(/NH:i/,@{$s});
    
    if(@grepres==1){
	my @t=split(":",$grepres[0]);
	$nh=$t[2]+0;
    }

    return $nh;
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts unique read ids from alignments.\n";
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
$parameters{"pathAlignments"}="NA";
$parameters{"forbiddenChromo"}="NA";
$parameters{"pathOutput"}="NA";
$parameters{"pathOutputCounts"}="NA";

my @defaultpars=("pathAlignments", "forbiddenChromo", "pathOutput", "pathOutputCounts");

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

my %forbiddenchr;
my @s=split(",", $parameters{"forbiddenChromo"});
foreach my $chr (@s){
    $forbiddenchr{$chr}=1;
}

my $nbf=keys %forbiddenchr;
print "Found ".$nbf." forbidden chromosomes: ".join(", ", keys %forbiddenchr)."\n";

my $path=$parameters{"pathAlignments"};

if(-e $path){
    print "Reading alignments from ".$path."\n";

    my $nbkept=0;

    open(my $output, ">".$parameters{"pathOutput"});

    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "sam"){
	open($input, $path);
    } else{
	if($ext eq "bam"){
	    open($input, "samtools view $path |");
	} else{
	    if($ext eq "gz"){
		open($input, "zcat $path |");
	    } else{
		print "Unknown file format\n";
		exit(1);
	    }
	}
    }

    my $line=<$input>;
  
    while($line){
	
	if($line eq ""){
	    $line=<$input>;
	    next;
	}
	
	my $firstchar=substr $line,0,1;
	
	if($firstchar eq "@"){
	    $line=<$input>;
	    next;
	}
	
	chomp $line;
	my @s=split("\t",$line);

	my $chr=$s[2];

	if(!exists $forbiddenchr{$chr}){
	
	    my $nbhits=findNbHits(\@s);
	    
	    if($nbhits==1){
		my $id=$s[0];
		print $output $id."\n";

		$nbkept++;
	    }
	}
	
	$line=<$input>;
    }

    close($output);
    close($input);


    open(my $outputcounts, ">".$parameters{"pathOutputCounts"});
    
    print $outputcounts $nbkept."\n";

    close($outputcounts);

}
else{
    print "Cannot find file: ".$path."\n";
}

print "Done.\n";

##############################################################
