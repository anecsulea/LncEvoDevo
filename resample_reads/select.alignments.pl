use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(uniq);
use strict;

##############################################################
##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script selects alignments.\n";
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
$parameters{"pathAllReads"}="NA";
$parameters{"pathSelectedReads"}="NA";
$parameters{"writeHeader"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathAlignments", "pathAllReads", "pathSelectedReads", "writeHeader", "pathOutput");

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

print "Reading all reads...\n";

my @allids;

my $pathall=$parameters{"pathAllReads"};
my @s=split("\\.",$pathall);
my $ext=$s[-1];

my $input;

if($ext eq "gz"){
    open($input, "zcat $pathall |");
} else{
    open($input, $pathall);
}

my $line=<$input>;

while($line){
    chomp $line;
    push(@allids, $line);
    
    $line=<$input>;
}

close($input);

my $nball=@allids;

print "Found ".$nball." reads in total.\n";

print "Done.\n";

##############################################################

print "Reading selected reads...\n";

my %selected;

open(my $input, $parameters{"pathSelectedReads"});

my $line=<$input>;

while($line){
    chomp $line;
    my $index=$line+0;
    my $id=$allids[$index];
    $selected{$id}=1;
    
    $line=<$input>;
}

close($input);

my $nbsel=keys %selected;

print "Found ".$nbsel." selected reads.\n";

print "Done.\n";

##############################################################

my $path=$parameters{"pathAlignments"};

if(-e $path){
    print "Reading alignments from ".$path."\n";

    open(my $output, ">".$parameters{"pathOutput"});

    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "sam"){
	open($input, $path);
    } else{
	if($ext eq "bam"){
	    open($input, "samtools view -h $path |");
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
	    if($parameters{"writeHeader"} eq "yes"){
		print $output $line; ## not chomped up yet
	    } else{
		print "not writing header lines.\n";
	    }
	    
	    $line=<$input>;
	    next;
	}
	
	chomp $line;
	my @s=split("\t",$line);

	my $id=$s[0];

	if(exists $selected{$id}){
            ## replace sequence and quality with "*", to make files lighter
	    $s[9]="*";
	    $s[10]="*"; 
	    print $output join("\t", @s)."\n";
	}

	$line=<$input>;
    }

    close($output);
    close($input);
}
else{
    print "Cannot find file: ".$path."\n";
}
print "Done.\n";

##############################################################
