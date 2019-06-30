use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#########################################################################

sub readAlignments{
    my $pathin=$_[0];
    my $species=$_[1]; ## required species
    my $minlen=$_[2]; ## minimum alignment length (excluding gaps, for each species)
    my $aln=$_[3];
    
    my $input;

    my @ext=split("\\.",$pathin);
    my $nbext=@ext;
    my $ext=$ext[$nbext-1];

    if($ext eq "gz"){
	open($input,"zcat $pathin | ");
    }
    else{
	open($input,$pathin);
    }

    my $line=<$input>;
    my $firstchar=substr $line,0,1;

    while($firstchar eq "#"){
	$line=<$input>;
	$firstchar=substr $line,0,1;
    }
 
    my $currentscore="NA";
    
    my $indexaln=0;

    while($line){

	chomp $line;
	$firstchar=substr $line,0,1;
	
	if($firstchar eq "a"){
	    
	    $indexaln++;

	    my @s=split(" ",$line);
	    my $score=$s[1];
	    my @t=split("=",$score);
	    $score=$t[1]+0;
	    $currentscore=$score;

	    $aln->{$indexaln}={};
	}

	if($firstchar eq "s"){
	    my @s=split(" ",$line);
	    
	    my $spchr=$s[1];
	    my @t=split("\\.",$spchr); 
	    my $sp=$t[0];
	    my $chr=$t[1];

	    my $start=$s[2]+0; ## 0-based
	    my $ungappedlen=$s[3]+0;
	    my $strand=$s[4];
	    my $chrlen=$s[5]+0;
	    my $sequence=$s[6];
	    
	    $aln->{$indexaln}{$sp}={"start"=>$start,"strand"=>$strand,"sequence"=>$sequence,"chrlen"=>$chrlen,"ungappedlen"=>$ungappedlen};
	    
	}
		
	$line=<$input>;
    }

    close($input);
 

    ## remove alignments that do not have all species and/or which are too short

    my @indexes=keys %{$aln};
    
    foreach my $index (@indexes){
	my $nbsp=keys %{$aln->{$index}};

	my $allin=1;
	my $alllong=1;

	foreach my $sp (@{$species}){
	    if(!(exists $aln->{$index}{$sp})){
		$allin=0;
		last;
	    }

	    my $len=$aln->{$index}{$sp}{"ungappedlen"};
	    
	    if($len<$minlen){
		$alllong=0;
		last;
	    }
	}

	if($allin==0 || $alllong==0){
	    delete $aln->{$index};
	}
    }
   
    my $nbkept=keys %{$aln};

    # print "Kept ".$nbkept." alignments.\n";
}

################################################################################

sub extractAlignmentStats{
    my $aln=$_[0];
    my $alnstats=$_[1];
  
    my @allindexes=keys %{$aln};
    my $firstindex=$allindexes[0];
    my @splist=keys %{$aln->{$firstindex}};
    my $firstsp=$splist[0];

    my $nbungapped=0;
    my $nbidentical=0;
    
    foreach my $index (keys %{$aln}){
	## now go over the alignment base by base

	my $alnlength=length $aln->{$index}{$firstsp}{"sequence"};

	for(my $i=0; $i<$alnlength; $i++){
	    my %bases;

	    my $allungap=1;
	    
	    foreach my $sp (keys %{$aln->{$index}}){
		my $base=uc (substr $aln->{$index}{$sp}{"sequence"}, $i, 1);

		if($base eq "-"){
		    $allungap=0;
		}
		else{
		    $bases{$base}=1;
		}
	    }
	    
	    if($allungap==1){
		$nbungapped++;

		my $nbdiffbases=keys %bases;

		if($nbdiffbases==1){
		    $nbidentical++;
		}
	    }
	}
    }

    $alnstats->{"ungapped"}=$nbungapped;
    $alnstats->{"identical"}=$nbidentical;
}

################################################################################

sub readClusters{
    my $pathin=$_[0];
    my $sp1=$_[1];
    my $sp2=$_[2];
    my $clusters=$_[3];
    
    open(my $input, $pathin);
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;
    
    my $index=0;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $g1=$s[$header{"Genes.".$sp1}];
	my $g2=$s[$header{"Genes.".$sp2}];
	
	if($g1 ne "NA" && $g2 ne "NA"){
	    my @genes1=split(";", $g1);
	    my @genes2=split(";", $g2);
	    
	    $clusters->{$index}={$sp1=>[], $sp2=>[]};

	    push(@{$clusters->{$index}{$sp1}}, @genes1);
	    push(@{$clusters->{$index}{$sp2}}, @genes2);
		
	    $index++;
	}
	
	$line=<$input>;
    }
        
    close($input);
}
################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts exon alignment stats from TBA alignments. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

################################################################################
################################################################################

my %parameters;
$parameters{"species1"}="NA";
$parameters{"species2"}="NA";
$parameters{"pathClusters"}="NA";
$parameters{"dirTBA"}="NA";
$parameters{"minAlignmentLength"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("species1", "species2", "pathClusters", "dirTBA", "minAlignmentLength", "pathOutput");

my @numericpars=();


my %numericpars;

foreach my $par (@numericpars){
    $numericpars{$par}=1;
}

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

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

#####################################################################
#####################################################################

my $sp1=$parameters{"species1"};
my $sp2=$parameters{"species2"};
my @species=($sp1, $sp2);

#####################################################################

print "Reading clusters...\n";

my %clusters;
readClusters($parameters{"pathClusters"}, $sp1, $sp2, \%clusters);

my $nbclust=keys %clusters;

print "Found ".$nbclust." non-empty clusters.\n";
print "Done.\n";

#####################################################################

print "Reading alignments and writing output...\n";

my $minlength=$parameters{"minAlignedLength"}+0;

print "minimum length: ".$minlength."\n";

open(my $output, ">".$parameters{"pathOutput"});
print $output "ID.".$sp1."\tID.".$sp2."\tLengthUngapped\tLengthIdentical\n";

foreach my $idclust (keys %clusters){
    foreach my $gene1 (@{$clusters{$idclust}{$sp1}}){
	foreach my $gene2 (@{$clusters{$idclust}{$sp2}}){
	    my $pathTBA=$parameters{"dirTBA"}."/".$gene1.".".$gene2.".maf.gz";
	    
	    if(-e $pathTBA){
		my %aln;
		readAlignments($pathTBA, \@species, $minlength, \%aln);
		
		my %alnstats;
		extractAlignmentStats(\%aln, \%alnstats);
		
		print $output $gene1."\t".$gene2."\t".$alnstats{"ungapped"}."\t".$alnstats{"identical"}."\n";
		
	    } else{
		print "Weird! cannot find ".$pathTBA." for ".$gene1." and ".$gene2."\n";
	    }
	}
    }
}

close($output);

print "Done.\n";

#####################################################################
#####################################################################
