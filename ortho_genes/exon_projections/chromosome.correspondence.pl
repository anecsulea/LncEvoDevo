use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
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

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts splice site coordinates from Ensembl\n";
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

$parameters{"pathEnsembl"}="NA";
$parameters{"pathUCSC"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathEnsembl","pathUCSC","pathOutput");
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

print "Reading chromosome sequences...\n";

my %ensembl;
readFasta($parameters{"pathEnsembl"}, \%ensembl);

my $nbens=keys %ensembl;

print "Found ".$nbens." chromosomes in Ensembl.\n";

my %ucsc;
readFasta($parameters{"pathUCSC"}, \%ucsc);

my $nbucsc=keys %ucsc;

print "Found ".$nbucsc." chromosomes in UCSC.\n";

print "Done.\n";

##############################################################

print "Extracting chromosome sizes...\n";

my %enssizes;

foreach my $chr (keys %ensembl){
    my $size=length $ensembl{$chr};

    if(exists $enssizes{$size}){
	push(@{$enssizes{$size}}, $chr);
    }
    else{
	$enssizes{$size}=[$chr];
    }
}

my %ucscsizes;

foreach my $chr (keys %ucsc){
    my $size=length $ucsc{$chr};

    if(exists $ucscsizes{$size}){
	push(@{$ucscsizes{$size}}, $chr);
    }
    else{
	$ucscsizes{$size}=[$chr];
    }
}

print "Done.\n";

##############################################################

print "Comparing chromosomes and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "Ensembl\tUCSC\tSize\n";

foreach my $size (keys %enssizes){
    if(exists $ucscsizes{$size}){
	
	foreach my $chrens (@{$enssizes{$size}}){
	    my %corresp;

	    my $seqens=uc $ensembl{$chrens};
	    
	    foreach my $chrucsc (@{$ucscsizes{$size}}){
		my $sequcsc=uc $ucsc{$chrucsc};

		if($seqens eq $sequcsc){
		    $corresp{$chrucsc}=1;
		}
	    }

	    my $nbcorr=keys %corresp;

	    if($nbcorr==1){
		my @corr=keys %corresp;

		print $output $chrens."\t".$corr[0]."\t".$size."\n";
	    }
	    else{
		print "Found ".$nbcorr."\t correspondences for ".$chrens."\n";

		if($nbcorr>1){
		    print join(",", keys %corresp)."\n";
		}
	    }
	}
    }
}

close($output);

print "Done.\n";

##############################################################
