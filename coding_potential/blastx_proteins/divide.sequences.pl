use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#############################################################################
#############################################################################

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

################################################################################

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

################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script divides large fasta files. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

#################################################################################
#################################################################################

my %parameters;

$parameters{"pathFastaInput"}="NA";
$parameters{"nbParts"}="NA";
$parameters{"prefixOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathFastaInput", "nbParts", "prefixOutput");

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

###################################################################################
###################################################################################

print "Reading sequences...\n";
my %sequences;
readFasta($parameters{"pathFastaInput"}, \%sequences);
print "Done.\n";

my $nbparts=$parameters{"nbParts"}+0;

print "Dividing sequences in ".$nbparts." files.\n";

my @ids=keys %sequences;
my $nbseq=@ids;


for(my $i=0; $i<$nbparts; $i++){
    open(my $output, ">".$parameters{"prefixOutput"}."_part".$i.".fa");

    for(my $k=0; $k<$nbseq; $k++){

	if(($k%$nbparts)==$i){
	    my $id=$ids[$k];
	    my $seq=$sequences{$id};

	    writeSequence($seq, $id, $output);
	}
    }

    close($output);
}

print "Done.\n";

###################################################################################
###################################################################################
