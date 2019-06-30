use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

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
	open($input,$path);
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

###############################################################

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

##############################################################

sub printHelp{
    
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts fake reads for each genome\n";
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
$parameters{"pathGenome"}="NA";
$parameters{"readLength"}="NA";
$parameters{"step"}="NA";
$parameters{"chr"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathGenome","readLength","step","chr","pathOutput");
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

print "Reading genome sequence...\n";
my %genome;
readFasta($parameters{"pathGenome"},\%genome);
my $nbchr=keys %genome;
print "There are ".$nbchr." sequences in ".$parameters{"pathGenome"}."\n";
print "Done.\n";


my $chr=$parameters{"chr"};

if(exists $genome{$chr}){

    print "Generating fake reads and writing output...\n";
    
    open(my $output, ">".$parameters{"pathOutput"});
    
    my $readlen=$parameters{"readLength"}+0;
    my $step=$parameters{"step"}+0;
    
    my $onlyN = "N" x $readlen;
    
    print "We don't output reads that look like this: ".$onlyN."\n";
    
    print "Read length ".$readlen." step ".$step."\n";
    
    my $index=0;
    
    if(exists $genome{$chr}){
	my $size=length $genome{$chr};
	
	print $chr." ".$size."\n";
	
	my $start=0;
	
	while($start<($size-$readlen)){
	    my $seq=substr $genome{$chr},$start,$readlen;
	    $seq=uc $seq;
	    
	    if($seq ne $onlyN){
		my $seqname="read".$index." ".$chr.":".$start."-".($start+$readlen-1);
		writeSequence($seq,$seqname,$output);
		$index++;
	    }
	    
	    $start+=$step;
	}
    }
    
    close($output);
    
    print "Done.\n";
    
}
else{
    print "Cannot find chromosome ".$chr.", not doing anything\n";
}

##############################################################
