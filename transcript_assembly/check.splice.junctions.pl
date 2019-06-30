use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $transcripts=$_[1];
    my $txgenes=$_[2];
    
    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon"){
	    my $chr=$s[0];
	    my $start=$s[3]+0;
	    my $end=$s[4]+0;
	    my $ss=$s[6];
	    my $strand="NA";
	    
	    if($ss eq "+"){
		$strand="1";
	    } else{
		if($ss eq "-"){
		    $strand="-1";
		} 
	    }

	   
	    ## we discard here unstranded transcripts
	    my $info=$s[8];
	    my @infoarray=split(";", $info);
	    
	    my $txid=findInfo("transcript_id", \@infoarray);
	    
	    if($txid eq "NA"){
		print "Weird! cannot find transcript info in ".$line."\n";
		exit(1);
	    }

	    my $gene=findInfo("gene_id", \@infoarray);

	    if($gene eq "NA"){
		print "Weird! cannot find gene info for ".$txid."\n";
		exit(1);
	    }

	    $txgenes->{$txid}=$gene;
	    
	    if(exists $transcripts->{$txid}){
		push(@{$transcripts->{$txid}{"start"}}, $start);
		push(@{$transcripts->{$txid}{"end"}}, $end);
	    } else{
		$transcripts->{$txid}={"start"=>[$start], "end"=>[$end], "chr"=>$chr, "strand"=>$strand};
	    }
	    
	}
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];
    
    my $res="NA";
    my $nbfound=0;
    
    my @grepres=grep(/${pattern}/,@{$array});
    
    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
	$nbfound++;
    }
    else{
	foreach my $g (@grepres){
	    my @u=split(" ",$g);
	    
	    if($u[0] eq $pattern){
		$nbfound++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }
    
    if($nbfound==1){
	return $res;
    } else{
	return "NA";
    }
}

##############################################################

sub readCorrectJunctions{
    my $pathin=$_[0];
    my $junctions=$_[1];
    
    open(my $input, $pathin);
    
    my $line=<$input>;
    chomp $line;
    my @s=split("\t",$line);
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[$header{"Chromosome"}];
	my $start=$s[$header{"Start"}];
	my $end=$s[$header{"End"}];
	my $strand=$s[$header{"ProbableStrand"}];
	
	my $key=$chr.",".$start.",".$end.",".$strand;
	
	$junctions->{$key}=1;
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub readWrongJunctions{
    my $pathin=$_[0];
    my $junctions=$_[1];
    
    open(my $input, $pathin);
    
    my $line=<$input>;
       
    while($line){
	chomp $line;
	my @s=split(" ", $line);
	my $chr=$s[0];
	my @t=split("_",$s[1]);
	my $start=$t[0];
	my $end=$t[1];
	my $strand=$t[2];

	my $splicestrand=$s[5];
	
	if($splicestrand ne "NA"){ ## otherwise could be just genome sequence errors
	    my $key=$chr.",".$start.",".$end.",".$strand;
	    $junctions->{$key}=1;
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
    print "This script checks splice junctions.\n";
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
$parameters{"pathAnnotGTF"}="NA";
$parameters{"pathCorrectJunctions"}="NA";
$parameters{"pathWrongJunctions"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathAnnotGTF", "pathCorrectJunctions", "pathWrongJunctions", "pathOutput");

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

print "Reading transcripts...\n";

my %transcripts;
my %txgenes;

readGTF($parameters{"pathAnnotGTF"}, \%transcripts, \%txgenes);

my $nbtx=keys %transcripts;
print "Found ".$nbtx." transcripts.\n";

print "Done.\n";

##############################################################

print "Reading junctions...\n";

my %correct;
my @pathscorrect=split(",", $parameters{"pathCorrectJunctions"});
foreach my $path (@pathscorrect){
    readCorrectJunctions($path, \%correct);
}

my %wrong;
my @pathswrong=split(",", $parameters{"pathWrongJunctions"});

foreach my $path (@pathswrong){
    readWrongJunctions($path, \%wrong);
}

my $nbcorrect=keys %correct;
my $nbwrong=keys %wrong;

print "Found ".$nbcorrect." correct junctions and ".$nbwrong." wrong junctions.\n";

print "Done.\n";

##############################################################

print "Analyzing transcripts and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "TranscriptID\tGeneID\tNbIntrons\tNbSupportedIntrons\tNbWrongStrand\tWrongJunctions\n";

foreach my $txid (keys %transcripts){
    my $gene=$txgenes{$txid};

    my $chr=$transcripts{$txid}{"chr"};
    my $strand=$transcripts{$txid}{"strand"};

    my $nbex=@{$transcripts{$txid}{"start"}};
      
    if($nbex>1){
	my $nbintrons=$nbex-1;
	my $nbcorrect=0;
	my $nbwrong=0;

	my %wr;

	for(my $i=0; $i<($nbex-1); $i++){
	    my $end1=${$transcripts{$txid}{"end"}}[$i];
	    my $start2=${$transcripts{$txid}{"start"}}[$i+1];

	    if($start2<=$end1){
		print "Weird! exons are not ordered!\n";
		exit(1);
	    }
	    
	    my $keyjunc=$chr.",".($end1+1).",".($start2-1).",".$strand;

	    if(exists $correct{$keyjunc}){
		$nbcorrect++;
	    }

	    if(exists $wrong{$keyjunc}){
		$nbwrong++;
		$wr{$keyjunc}=1;
	    }
	}
	
	if($nbwrong>0){
	    print $output $txid."\t".$gene."\t".$nbintrons."\t".$nbcorrect."\t".$nbwrong."\t".join(";",keys %wr)."\n";
	}
	else{
	    print $output $txid."\t".$gene."\t".$nbintrons."\t".$nbcorrect."\t".$nbwrong."\tNA\n";
	}
    } else{
	print $output $txid."\t".$gene."\t0\t0\t0\tNA\n";
    }
}

close($output);

print "Done.\n";

##############################################################
