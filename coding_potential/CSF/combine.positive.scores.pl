use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX;
use strict;

##############################################################

sub readPositiveScores{

    my $pathin=$_[0];
    my $scores=$_[1];

    open(my $input, $pathin);
    
    my %header;
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	my $sumscore=$s[$header{"SumScore"}]+0.0;
	my $nbinfo=$s[$header{"NbInfoCodons"}]+0;

	if(!(exists $scores->{$chr})){
	    $scores->{$chr}={};
	}
	
	if(!(exists $scores->{$chr}{$strand})){
	    $scores->{$chr}{$strand}={"start"=>[], "end"=>[], "sumscore"=>[], "nbinfo"=>[]};
	}
	
	push(@{$scores->{$chr}{$strand}{"start"}}, $start);
	push(@{$scores->{$chr}{$strand}{"end"}}, $end);
	push(@{$scores->{$chr}{$strand}{"sumscore"}}, $sumscore);
	push(@{$scores->{$chr}{$strand}{"nbinfo"}}, $nbinfo); 
       
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub orderScores{
    my $unordered=$_[0];
    my $ordered=$_[1];

    foreach my $chr (keys %{$unordered}){
	$ordered->{$chr}={};

	foreach my $strand (keys %{$unordered->{$chr}}){
	    $ordered->{$chr}{$strand}={"start"=>[], "end"=>[], "sumscore"=>[], "nbinfo"=>[]};

	    my %hashpos;
	    my $nbstart=@{$unordered->{$chr}{$strand}{"start"}};

	    for(my $i=0; $i<$nbstart; $i++){
		my $start=${$unordered->{$chr}{$strand}{"start"}}[$i];
		my $end=${$unordered->{$chr}{$strand}{"end"}}[$i];
		my $sumscore=${$unordered->{$chr}{$strand}{"sumscore"}}[$i];
		my $nbinfo=${$unordered->{$chr}{$strand}{"nbinfo"}}[$i];

		if(exists $hashpos{$start}){
		    if(exists $hashpos{$start}{$end}){
			push(@{$hashpos{$start}{$end}{"sumscore"}}, $sumscore);
			push(@{$hashpos{$start}{$end}{"nbinfo"}}, $nbinfo);
		    }
		    else{
			$hashpos{$start}{$end}={"sumscore"=>[$sumscore], "nbinfo"=>[$nbinfo]};
		    }
		}
		else{
		    $hashpos{$start}={$end=>{"sumscore"=>[$sumscore], "nbinfo"=>[$nbinfo]}};
		}
	    }

	    my @uniquestart=keys %hashpos;
	    my @sortedstart=sort {$a<=>$b} @uniquestart;

	    foreach my $start (@sortedstart){
		my @uniqueend=keys %{$hashpos{$start}};
		my @sortedend=sort {$a<=>$b} @uniqueend;

		foreach my $end (@sortedend){
		    my $nb=@{$hashpos{$start}{$end}{"sumscore"}};

		    for(my $i=0; $i<$nb; $i++){
			my $sumscore=${$hashpos{$start}{$end}{"sumscore"}}[$i];
			my $nbinfo=${$hashpos{$start}{$end}{"nbinfo"}}[$i];
			
			push(@{$ordered->{$chr}{$strand}{"start"}}, $start);
			push(@{$ordered->{$chr}{$strand}{"end"}}, $end);
			push(@{$ordered->{$chr}{$strand}{"sumscore"}}, $sumscore);
			push(@{$ordered->{$chr}{$strand}{"nbinfo"}}, $nbinfo);
		    }
		}
	    }
	}
    }
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script combines positive CSF scores.\n";
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

$parameters{"pathsScores"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathsScores", "pathOutput");


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

print "Reading positive CSF scores...\n";
my %scores;

my @paths=split(",", $parameters{"pathsScores"});

foreach my $path (@paths){
    print $path."\n";
 
    readPositiveScores($path, \%scores);
}

print "Done.\n";

##############################################################

print "Ordering scores...\n";

my %orderedscores;

orderScores(\%scores, \%orderedscores);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "Chr\tStart\tEnd\tStrand\tSumScore\tNbInfoCodons\n";

foreach my $chr (keys %orderedscores){
    foreach my $strand (keys %{$orderedscores{$chr}}){
	my $nb=@{$orderedscores{$chr}{$strand}{"sumscore"}};
	
	for(my $i=0; $i<$nb; $i++){
	    my $start=${$orderedscores{$chr}{$strand}{"start"}}[$i];
	    my $end=${$orderedscores{$chr}{$strand}{"end"}}[$i];
	    my $sumscore=${$orderedscores{$chr}{$strand}{"sumscore"}}[$i];
	    my $nbinfo=${$orderedscores{$chr}{$strand}{"nbinfo"}}[$i];

	    print $output $chr."\t".$start."\t".$end."\t".$strand."\t".$sumscore."\t".$nbinfo."\n";
	}
    }
}

close($output);

print "Done.\n";

##############################################################
