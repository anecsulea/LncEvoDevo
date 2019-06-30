use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readProjectedExons{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $genesexons=$_[2];
    
    open(my $input, $pathin);
    
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $geneid=$s[$header{"GeneID"}];
		
	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0; ## 1-based
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	
	my $exonid=$chr.",".$start.",".$end.",".$strand;

	$exoncoords->{$exonid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	
	if(exists $genesexons->{$geneid}){
	    $genesexons->{$geneid}{$exonid}=1;
	}
	else{
	    $genesexons->{$geneid}={$exonid=>1};
	}
        
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub makeBlocks{
    my $genesexons=$_[0];
    my $exoncoords=$_[1];
    my $genesblocks=$_[2];
    my $blockcoords=$_[3];

    foreach my $gene (keys %{$genesexons}){
	my %hashpos;
	
	my $chr="NA";
	my $strand="NA";

	$genesblocks->{$gene}={};

	foreach my $exonid (keys %{$genesexons->{$gene}}){
	    my $start=$exoncoords->{$exonid}{"start"};
	    my $end=$exoncoords->{$exonid}{"end"};
	    
	    if($chr eq "NA"){
		$chr=$exoncoords->{$exonid}{"chr"};
		$strand=$exoncoords->{$exonid}{"strand"};
	    } else{
		if($exoncoords->{$exonid}{"chr"} ne $chr){
		    print "Saw two different chromosomes for ".$gene.": ".$chr." ".$exoncoords->{$exonid}{"chr"}."\n";
		    exit(1);
		}

		if($exoncoords->{$exonid}{"strand"} ne $strand){
		    print "Saw two different strands for ".$gene.": ".$strand." ".$exoncoords->{$exonid}{"strand"}."\n";
		    exit(1);
		}
	    }
	    
	    if(exists $hashpos{$start}){
		if($end>$hashpos{$start}){
		    $hashpos{$start}=$end;
		}
	    }
	    else{
		$hashpos{$start}=$end;
	    }
	}
    

	my @uniquepos=keys %hashpos;
	my @sortedpos=sort {$a<=>$b} @uniquepos;
	
	my $currentstart=$sortedpos[0];
	my $currentend=$hashpos{$sortedpos[0]};
	
	for(my $i=1; $i<@sortedpos; $i++){
	    my $thisstart=$sortedpos[$i];
	    my $thisend=$hashpos{$sortedpos[$i]};
	    
	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		if($currentend<$thisend){
		    $currentend=$thisend;
		}
	    }
	    else{
		my $currentid=$chr.",".$currentstart.",".$currentend.",".$strand;
		$blockcoords->{$currentid}={"chr"=>$chr, "start"=>$currentstart, "end"=>$currentend, "strand"=>$strand};
		$genesblocks->{$gene}{$currentid}=1;
		
		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}
	
	my $currentid=$chr.",".$currentstart.",".$currentend.",".$strand;
	$blockcoords->{$currentid}={"chr"=>$chr, "start"=>$currentstart, "end"=>$currentend, "strand"=>$strand};
	$genesblocks->{$gene}{$currentid}=1;
	
    }
}

##############################################################

sub writeGTF{
    my $genesblocks=$_[0];
    my $blockcoords=$_[1];
    my $pathout=$_[2];

    open(my $output, ">".$pathout);

    foreach my $geneid (keys %{$genesblocks}){
	foreach my $exonid (keys %{$genesblocks->{$geneid}}){
	    my $chr=$blockcoords->{$exonid}{"chr"};
	    my $start=$blockcoords->{$exonid}{"start"};
	    my $end=$blockcoords->{$exonid}{"end"};
	    my $strand=$blockcoords->{$exonid}{"strand"};

	    my $outstrand="NA";

	    if($strand eq "1"){
		$outstrand="+";
	    } else{
		if($strand eq "-1"){
		    $outstrand="-";
		} else{
		    print "Weird strand for ".$geneid." and ".$exonid.": ".$strand."\n";
		}
	    }

	    my $features="gene_id \"".$geneid."\"; transcript_id \"".$geneid.".tx\";";

	    print $output $chr."\tProjection\texon\t".$start."\t".$end."\t.\t".$outstrand."\t.\t".$features."\n";
	}
    }

    close($output);
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script writes projected exon coordinates in GTF format\n";
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
$parameters{"pathProjectedExons"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathProjectedExons", "pathOutput");

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

print "Reading projected annotations...\n";

my %exoncoords;
my %genesexons;

readProjectedExons($parameters{"pathProjectedExons"}, \%exoncoords, \%genesexons);

my $nbex=keys %exoncoords;
my $nbg=keys %genesexons;

print "Found ".$nbex." exons and ".$nbg." genes.\n";

print "Done.\n";

##############################################################

print "Making exon blocks...\n";

my %blockcoords;
my %genesblocks;

makeBlocks(\%genesexons, \%exoncoords, \%genesblocks, \%blockcoords);
  
print "Done.\n";

##############################################################

print "Writing GTF output...\n";

writeGTF(\%genesblocks, \%blockcoords, $parameters{"pathOutput"});

print "Done.\n";

##############################################################
