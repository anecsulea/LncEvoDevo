use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readExonCoords{
    my $pathExons=$_[0];
    my $exons=$_[1];

    open(my $input,$pathExons);
    
    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);
    
    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $origstart=$s[$header{"NonoverlapStart"}];
	my $origend=$s[$header{"NonoverlapEnd"}];

	if($origstart ne "NA" && $origend ne "NA"){
	    my $geneid=$s[$header{"GeneID"}];
	    my $chr=$s[$header{"Chr"}];
	    
	    my $start=$s[$header{"NonoverlapStart"}]+0;
	    my $end=$s[$header{"NonoverlapEnd"}]+0;
	    my $strand=$s[$header{"Strand"}];
	    
	    if(exists $exons->{$geneid}){
		push(@{$exons->{$geneid}{"start"}}, $start);
		push(@{$exons->{$geneid}{"end"}}, $end);
	    }
	    else{
		$exons->{$geneid}={"chr"=>$chr,"strand"=>$strand, "start"=>[$start], "end"=>[$end]};
	    }
	}

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub makeExonBlocks{
    my $exoncoords=$_[0];
    my $collapse=$_[1];
    my $exonblocks=$_[2];

    foreach my $gene (keys %{$exoncoords}){

	my $chr=$exoncoords->{$gene}{"chr"};
	my $strand=$exoncoords->{$gene}{"strand"};

	my %hashcoords;

	my $nbexons=@{$exoncoords->{$gene}{"start"}};

	for(my $i=0; $i<$nbexons; $i++){
	    my $b=${$exoncoords->{$gene}{"start"}}[$i];
	    my $e=${$exoncoords->{$gene}{"end"}}[$i];
	    
	    if(exists $hashcoords{$b}){
		if($e>$hashcoords{$b}){
		    $hashcoords{$b}=$e;
		}
	    }
	    else{
		$hashcoords{$b}=$e;
	    }
	}
	
	$exonblocks->{$gene}={"chr"=>$chr,"strand"=>$strand,"start"=>[],"end"=>[]};
	
	my @uniquestart = keys %hashcoords;
	my @sortedstart = sort {$a <=> $b} @uniquestart;

	my $nbex=@sortedstart;
	
	my $currentstart=$sortedstart[0];
	my $currentend=$hashcoords{$currentstart};

	for(my $u=1; $u<$nbex; $u++){
	    my $thisstart=$sortedstart[$u];
	    my $thisend=$hashcoords{$thisstart};
	    
	    ## cluster blocks if they overlap
	    
	    if($thisstart>=$currentstart && $thisstart<=($currentend+$collapse+1)){  
		
		## we only change the end if it's larger than the current position
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    }
	    else{
		push(@{$exonblocks->{$gene}{"start"}},$currentstart);
		push(@{$exonblocks->{$gene}{"end"}},$currentend);
			
		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	    
	}
	
	## don't forget the last block
	
	push(@{$exonblocks->{$gene}{"start"}},$currentstart);
	push(@{$exonblocks->{$gene}{"end"}},$currentend);
	
    }
}

##############################################################

sub writeExons{

    my $refblocks=$_[0];
    my $pathoutput=$_[1];

    open(my $output,">".$pathoutput);
    
    foreach my $gene (keys %{$refblocks}){
	my $chr=$refblocks->{$gene}{"chr"};
	my $strand=$refblocks->{$gene}{"strand"};
	
	my $nbblocks=@{$refblocks->{$gene}{"start"}};
	
	for(my $i=0;$i<$nbblocks;$i++){
	    my $exonid=$gene.".".($i+1);
	    
	    my $start=${$refblocks->{$gene}{"start"}}[$i];
	    my $end=${$refblocks->{$gene}{"end"}}[$i];
	  
	    
	    print $output $gene."\t".$exonid."\t".$chr."\t".$start."\t".$end."\t".$strand."\n";
	}
	
    }
        
    close($output);
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script constructs exon blocks from non-overlapping exon coords.\n";
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

$parameters{"pathExonCoords"}="NA";
$parameters{"collapseDistance"}=10;
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonCoords","collapseDistance","pathOutput");


my %numericpars;
my @numericpars=("collapseDistance");

foreach my $par (@numericpars){
    $numericpars{$par}=1;
}

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

#####################################################################
#####################################################################

print "Reading exon coordinates...\n";
my %exoncoords;
readExonCoords($parameters{"pathExonCoords"},\%exoncoords);
print "Done.\n\n";

my $nbg=keys %exoncoords;
print "There are ".$nbg." genes.\n";

#####################################################################

print "Making exon blocks...\n";
my %exonblocks;
my $collapse=$parameters{"collapseDistance"}+0;

print "Collapsing blocks that are less than ".$collapse."bp apart.\n";

makeExonBlocks(\%exoncoords, $collapse, \%exonblocks);
print "Done.\n\n";

#####################################################################

print "Writing output...\n";
writeExons(\%exonblocks,$parameters{"pathOutput"});
print "Done.\n\n";

#####################################################################
