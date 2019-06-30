use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readBED{
    my $pathin=$_[0];
    my $coords=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[0];
	my $prefix=substr $chr, 0,3;

	if($prefix eq "chr"){ ## otherwise we don't analyze it
	    $chr=substr $chr, 3;
	    my $start=$s[1]+1; ## 1-based, included
	    my $end=$s[2]+0; ## 1-based, included
	    my $id=$s[3];
	    my $strand=$s[5];
	    
	    if(exists $coords->{$id}){
		print "Weird! already saw ".$id."\n";
		exit(1);
	    } else{
		$coords->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	    }
	    
	}
	
	$line=<$input>;
    }

    close($input);
}

###################################################################

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

############################################################################

sub reverseComplement{
    my $sequence=$_[0];
    
    my $rev=reverse $sequence;

    $rev=~s/A/X/g;
    $rev=~s/C/Y/g;
    $rev=~s/G/Z/g;
    $rev=~s/T/W/g;

    $rev=~s/X/T/g;
    $rev=~s/Y/G/g;
    $rev=~s/Z/C/g;
    $rev=~s/W/A/g;

    return $rev;
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

###############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts lifted promoter sequences. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

############################################################################
############################################################################

my %parameters;

$parameters{"pathCoordinates"}="NA";
$parameters{"pathGenomeSequence"}="NA";
$parameters{"outprefix"}="NA";
$parameters{"pathOutputSequences"}="NA";

my @defaultpars=("pathCoordinates", "pathGenomeSequence", "outprefix", "pathOutputSequences");

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

##############################################################################
##############################################################################

print "Reading lifted promoter coordinates...\n";

my %coords;
readBED($parameters{"pathCoordinates"}, \%coords);

my $nblift=keys %coords;
print "Found ".$nblift." coordinates.\n";

print "Done.\n";

##############################################################################

print "Reading genome sequence...\n";

my %genome;
readFasta($parameters{"pathGenomeSequence"}, \%genome);

my $nbchr=keys %genome;

print "Saw ".$nbchr." chromosomes.\n";

print "Done.\n";

##############################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutputSequences"});
my $prefix=$parameters{"outprefix"};

print "Out prefix: ".$prefix."\n";

foreach my $id (keys %coords){
    my $chr=$coords{$id}{"chr"};

    if(exists $genome{$chr}){
	my $start=$coords{$id}{"start"};
	my $end=$coords{$id}{"end"};
	my $strand=$coords{$id}{"strand"};

	my $len=$end-$start+1;

	if($len>=10){
	    
	    my $seq=substr $genome{$chr}, ($start-1), ($end-$start+1);
	    
	    if($strand eq "-"){
		$seq=reverseComplement($seq);
	    }
	    
	    my $outid=$prefix."_".$id;
	    
	    writeSequence($seq, $outid, $output);
	    
	}
    }
}

close($output);

print "Done.\n";


##############################################################################
