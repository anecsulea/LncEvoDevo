use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX;
use strict;

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script divides MAF alignments by chromosome.\n";
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

$parameters{"pathAlignment"}="NA";
$parameters{"refSpecies"}="NA";
$parameters{"chromosomes"}="NA";
$parameters{"dirOutput"}="NA";

my @defaultpars=("pathAlignment", "refSpecies", "chromosomes", "dirOutput");


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

print "Extracting accepted chromosomes...\n";

my @s=split(",", $parameters{"chromosomes"});

my %chromo;

foreach my $chr (@s){
    $chromo{$chr}=1;
}

my $nbchr=keys %chromo;

print "There are ".$nbchr." accepted chromosomes.\n";

print "Done.\n";

##############################################################

print "Preparing output files...\n";

my %outputfiles;

foreach my $chr (keys %chromo){
    my $pathout=$parameters{"dirOutput"}.$chr.".maf";
    my $output;
    open($output, ">".$pathout);
    $outputfiles{$chr}=$output;
}

print "Done.\n";

##############################################################

print "Reading alignment and writing output...\n";

my $pathin=$parameters{"pathAlignment"};
my $input;

my @s=split("\\.", $pathin);
my $ext=$s[-1];

if($ext eq "gz"){
    open($input, "zcat $pathin |");
} else{
    open($input, $pathin);
}

my $line=<$input>;

my @currentaln;
my $currentchr="NA";

my $refsp=$parameters{"refSpecies"};

print "reference species: ".$refsp."\n";

my $nbread=0;

while($line){
    chomp $line;
    my $firstchar=substr $line,0,1;
    
    if($firstchar eq "a"){
	## new alignment

	if($currentchr ne "NA"){
	    ## if the reference species is not present in the alignment, we don't write output

	    if(exists $chromo{$currentchr}){
		my $output=$outputfiles{$currentchr};
		
		foreach my $line (@currentaln){
		    print $output $line."\n";
		}
	    }
	    
	    undef @currentaln;
	}

	$currentchr="NA"; ## for now we don't know the chromosome
    }
    
    if($firstchar eq "s"){
	my @t=split(" ",$line);
	my $spchr=$t[1];
	my @u=split("\\.",$spchr);
	
	if($u[0] eq $refsp){
	    $currentchr=$u[1];
	}
    }
    
    push (@currentaln, $line);
    
    $line=<$input>;

    $nbread++;

    if($nbread%100000==0){
	print "read ".$nbread." lines.\n";
    }
}

## don't forget last alignment

if($currentchr ne "NA"){
    if(exists $chromo{$currentchr}){
	my $output=$outputfiles{$currentchr};
	
	foreach my $line (@currentaln){
	    print $output $line."\n";
	}
    }
    
    undef @currentaln;
}

close($input);

print "Done.\n";

##############################################################

print "Closing output files...\n";

foreach my $chr (keys %chromo){
    my $output=$outputfiles{$chr};
    close($output);
}

print "Done.\n";

##############################################################
