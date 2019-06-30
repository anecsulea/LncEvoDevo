use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readGeneticCode{
    my $pathin=$_[0];
    my $code=$_[1];

    open(my $input,$pathin);
    my $line=<$input>;

    while($line){
	my @s=split(" ",$line);
	$code->{$s[0]}=$s[1];

	$line=<$input>;
    }

    close($input);
}


##############################################################

sub readMatrices{
    my $path=$_[0];
    my $matrices=$_[1];

    open(my $input, $path);

    my $line=<$input>;

    my $sp1="NA";
    my $sp2="NA";
    
    while($line){
	chomp $line;

	my $firstchar=substr $line, 0, 1; 

	if($firstchar eq "#"){
	    my $species=substr $line,2;
	    my @s=split(" ",$species);
	    $sp1=$s[0];
	    $sp2=$s[1];
	}
	else{
	    if($sp1 eq "NA" || $sp2 eq "NA"){
		print "Weird!! NA species at ".$line."\n";
		exit(1);
	    }
	    
	    my @s=split("\t",$line);
	    my $codon1=$s[0];
	    my $codon2=$s[2];
	    my $nb=$s[4]+0;
	    
	    if(exists $matrices->{$sp1}){
		if(exists $matrices->{$sp1}{$sp2}){
		    if(exists $matrices->{$sp1}{$sp2}{$codon1}){
			if(exists $matrices->{$sp1}{$sp2}{$codon1}{$codon2}){
			    $matrices->{$sp1}{$sp2}{$codon1}{$codon2}+=$nb;
			}
			else{
			    $matrices->{$sp1}{$sp2}{$codon1}{$codon2}=$nb;
			}
		    }
		    else{
			$matrices->{$sp1}{$sp2}{$codon1}={$codon2=>$nb};
		    }
		}
		else{
		     $matrices->{$sp1}{$sp2}={$codon1=>{$codon2=>$nb}};
		 }
	    }
	    else{
		$matrices->{$sp1}={$sp2=>{$codon1=>{$codon2=>$nb}}};
	    }
	}

	$line=<$input>;

    }

    close($input);
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script combines CSF matrices.\n";
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
$parameters{"pathsMatrices"}="NA";
$parameters{"pathGeneticCode"}="NA";
$parameters{"pathOutput"}="NA";


my @defaultpars=("pathsMatrices","pathGeneticCode","pathOutput");


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

print "Reading matrices...\n";

my $paths=$parameters{"pathsMatrices"};
my $lastchar=substr $paths, -1;

if($lastchar eq ","){
    chop $paths;
}

my @paths=split(",", $paths);

my %matrices;

foreach my $path (@paths){

    if(-e $path){
	print $path."\n";
	readMatrices($path, \%matrices);
    }
    else{
	print "cannot find path ".$path."\n";
    }
}

print "Done.\n";
##############################################################

print "Reading genetic code...\n";

my %geneticcode;
readGeneticCode($parameters{"pathGeneticCode"},\%geneticcode);
my $nbcodons=keys %geneticcode;
print "There are ".$nbcodons." codons in the genetic code.\n";

print "Done.\n\n";

##############################################################

my @codons=keys %geneticcode;
    
print "Writing output...\n";

open(my $output,">".$parameters{"pathOutput"});

foreach my $sp1 (keys %matrices){
    foreach my $sp2 (keys %{$matrices{$sp1}}){
	print $output "# ".$sp1." ".$sp2."\n";
	foreach my $codon1 (@codons){
	    foreach my $codon2 (@codons){
		if(exists $matrices{$sp1}{$sp2}{$codon1}{$codon2}){
		    print $output $codon1."\t".$geneticcode{$codon1}."\t".$codon2."\t".$geneticcode{$codon2}."\t".$matrices{$sp1}{$sp2}{$codon1}{$codon2}."\n";
		}
		else{
		    print $output $codon1."\t".$geneticcode{$codon1}."\t".$codon2."\t".$geneticcode{$codon2}."\t"."\t0\n";
		}
	    }
	}
    }
}

close($output);

print "Done.\n\n";

##############################################################
