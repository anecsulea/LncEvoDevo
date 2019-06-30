use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $genes=$_[1];
    
    open(my $input, $pathin);
    
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $geneid=$s[0];
		
	my $chr=$s[2];
	my $start=$s[3]+0; ## 1-based
	my $end=$s[4]+0;
	my $strand=$s[5];
	
	
	my $exonid=$chr.",".$start.",".$end.",".$strand;
	
	if(exists $genes->{$geneid}){
	    $genes->{$geneid}{$exonid}=1;
	}
	else{
	    $genes->{$geneid}={$exonid=>1};
	}
        
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readProjectedExons{
    my $pathin=$_[0];
    my $exons=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;

    my %multiple;

    while($line){
	chomp $line;
	
	my @s=split("\t", $line);
	my $chr=$s[0];
	my $start=$s[1]+1; ## 1-based coordinates, included
	my $end=$s[2]+0; ## 1-based coordinate, included
	my $id=$s[3];
	my $strand=$s[5];
	
	if(exists $exons->{$id}){
	    $multiple{$id}=1;
	} 
	else{
	    $exons->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	}
	
	$line=<$input>;
    }
    
    close($input);
	
    my $nbmulti=keys %multiple;

    print "Found ".$nbmulti." exons with multiple (or split) mapping.\n";

    foreach my $id (keys %multiple){
	delete $exons->{$id};
    }
}

##############################################################

sub filterProjectedExons{
    my $projexons=$_[0];
    my $minsizeratio=$_[1];
    my $maxsizeratio=$_[2];
    my $filteredexons=$_[3];
    my $outLog=$_[4];
    
     my $nbbadsize=0;

    foreach my $id (keys %{$projexons}){
	my @s=split(",", $id);
	my $chr=$s[0];
	my $start=$s[1];
	my $end=$s[2];
	my $strand=$s[3];
	
	my $newchr=$projexons->{$id}{"chr"};
	my $newstart=$projexons->{$id}{"start"};
	my $newend=$projexons->{$id}{"end"};
	my $newstrand=$projexons->{$id}{"strand"};

	
	## we initialize projections at the expected positions 
	## if transcript boundaries, we don't check further
		
	my $oldsize=$end-$start+1.0;
	my $newsize=$newend-$newstart+1.0;
	
	my $sizeratio=$newsize/$oldsize;
	
	if($sizeratio>=$minsizeratio && $sizeratio<=$maxsizeratio){
	    $filteredexons->{$id}={"chr"=>$newchr, "start"=>$newstart, "end"=>$newend, "strand"=>$newstrand};
	}
	else{
	    $nbbadsize++;
	    print $outLog $id." bad_size ".$oldsize." ".$newsize." ".$sizeratio."\n";
	}
	
    }

    print "Rejected ".$nbbadsize." exons with bad size ratios.\n";
    
}

##############################################################

sub readChromosomeCorrespondence{
    my $pathin=$_[0]; 
    my $corresp12=$_[1]; 
    my $corresp21=$_[2]; 

    open(my $input, $pathin);

    my $line=<$input>; ## header
    
    my %dupliens;
    my %dupliucsc;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $ens=$s[0];
	my $ucsc=$s[1];
	
	if(exists $corresp12->{$ens}){
	    $dupliens{$ens}=1;
	    $dupliucsc{$ucsc}=1;
	}
	else{
	    $corresp12->{$ens}=$ucsc;
	}
	
	if(exists $corresp21->{$ucsc}){
	    $dupliens{$ens}=1;
	    $dupliucsc{$ucsc}=1;
	}
	else{
	    $corresp21->{$ucsc}=$ens;
	}
	
	$line=<$input>;
    }
    
    close($input);
    
    my $nbdupliens=keys %dupliens;
    
    if($nbdupliens>0){
	print "Found ".$nbdupliens." ambiguous Ensembl chromosomes: ".join("; ",keys %dupliens).".\n";

	foreach my $idens (keys %dupliens){
	    delete $corresp12->{$idens};
	}
    }

    my $nbdupliucsc=keys %dupliucsc;
    
    if($nbdupliucsc>0){
	print "Found ".$nbdupliucsc." ambiguous UCSC chromosomes: ".join("; ", keys %dupliucsc).".\n";

	foreach my $iducsc (keys %dupliucsc){
	    delete $corresp21->{$iducsc};
	}
    }
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script filters projected annotations.\n";
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

$parameters{"pathExonBlocks"}="NA";
$parameters{"pathProjectedExons"}="NA";
$parameters{"pathChromosomeCorrespondence"}="NA";
$parameters{"minSizeRatio"}="NA";
$parameters{"maxSizeRatio"}="NA";
$parameters{"pathOutputFilteredExons"}="NA";
$parameters{"pathOutputRejectedExons"}="NA";
$parameters{"pathOutputLog"}="NA";

my @defaultpars=("pathExonBlocks", "pathProjectedExons", "pathChromosomeCorrespondence", "minSizeRatio", "maxSizeRatio", "pathOutputFilteredExons", "pathOutputRejectedExons", "pathOutputLog");

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

print "Reading chromosome correspondence for reference species...\n";

my %chrcorresp12;
my %chrcorresp21;
readChromosomeCorrespondence($parameters{"pathChromosomeCorrespondence"}, \%chrcorresp12, \%chrcorresp21);

print "Done.\n";

##############################################################

print "Reading GTF...\n";

my %genes;

readExonBlocks($parameters{"pathExonBlocks"},  \%genes);

print "Done.\n";


##############################################################

print "Reading projected exons...\n";

my %projectedexons;

readProjectedExons($parameters{"pathProjectedExons"}, \%projectedexons);

my $nbproj=keys %projectedexons;

print "Found ".$nbproj." projected exons.\n";

print "Done.\n";

##############################################################

print "Filtering projected exons...\n";

my $minratio=$parameters{"minSizeRatio"}+0.0;
my $maxratio=$parameters{"maxSizeRatio"}+0.0;

print "Size ratio: ".$minratio." ".$maxratio."\n";

open(my $outLog, ">".$parameters{"pathOutputLog"});

my %filteredexons;

filterProjectedExons(\%projectedexons,  $minratio, $maxratio, \%filteredexons, $outLog);
  
my $nbkept=keys %filteredexons;

print "Kept ".$nbkept." exons.\n";

close($outLog);

print "Done.\n";

##############################################################

print "Writing output for filtered exons...\n";

open(my $output, ">".$parameters{"pathOutputFilteredExons"});

print $output "GeneID\tExonID\tChr\tStart\tEnd\tStrand\n";

foreach my $gene (keys %genes){
    foreach my $id (keys %{$genes{$gene}}){
	if(exists $filteredexons{$id}){
	    my $newchr=$filteredexons{$id}{"chr"};

	    if(exists $chrcorresp21{$newchr}){
		my $outchr=$chrcorresp21{$newchr};
		

		my $newstrand="NA";
		
		if($filteredexons{$id}{"strand"} eq "+"){
		    $newstrand="1";
		}
		else{
		    if($filteredexons{$id}{"strand"} eq "-"){
			$newstrand="-1";
		    } else{
			if($filteredexons{$id}{"strand"} eq "-1" || $filteredexons{$id}{"strand"} eq "1"){
			    $newstrand=$filteredexons{$id}{"strand"};
			} else{
			    print "Weird strand for ".$id." ".$filteredexons{$id}{"strand"} ."!!!\n";
			    exit(1);
			}
		    }
		}
		
		print $output $gene."\t".$id."\t".$outchr."\t".$filteredexons{$id}{"start"}."\t".$filteredexons{$id}{"end"}."\t".$newstrand."\n";
	    }
	    else{
		print "Not writing output for ".$id.", ambiguous chromosome: UCSC ".$newchr.".\n";
	    }
	}
    }
}
close($output);

print "Done.\n";

##############################################################

print "Writing output for rejected exons...\n";

open(my $output, ">".$parameters{"pathOutputRejectedExons"});

print $output "GeneID\tExonID\tChr\tStart\tEnd\tStrand\n";

foreach my $gene (keys %genes){
    foreach my $id (keys %{$genes{$gene}}){
	if((exists $projectedexons{$id}) && (!exists $filteredexons{$id})){
	    print $output $gene."\t".$id."\t".$projectedexons{$id}{"chr"}."\t".$projectedexons{$id}{"start"}."\t".$projectedexons{$id}{"end"}."\t".$projectedexons{$id}{"strand"}."\n";
	}
    }
}
close($output);

print "Done.\n";

##############################################################
