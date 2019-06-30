use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readBlocks{

    my $path=$_[0];
    my $refblocks=$_[1];
    ## format: tab-separated, gene, exon,  chr, strand, begin, end, tag
    
    open(my $input,$path);

    my $line=<$input>;

    while($line){
	
	chomp $line;
	my @s=split("\t",$line);
	my $gene=$s[0];
	my $chr=$s[2];
	my $begin=$s[3]+0;
	my $end=$s[4]+0;
	my $strand=$s[5];
	
	if(exists $refblocks->{$gene}){
	    push(@{$refblocks->{$gene}{"begin"}},$begin);
	    push(@{$refblocks->{$gene}{"end"}},$end);
	}
	else{
	    $refblocks->{$gene}={"chr"=>$chr,"strand"=>$strand,"begin"=>[$begin], "end"=>[$end]};
	}

	$line=<$input>;
    }

    close($input);
    
}

##############################################################

sub readCoverage{
    my $pathin=$_[0];
    my $cov=$_[1];

    my @s=split("\\.",$pathin);
    my $ext=$s[-1];
    
    my $input;
    
    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    }
    else{
	open($input, $pathin);
    }
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	my $val=$s[3]+0;

	if($val>0){
	    my $chr=$s[0];
	    my $start=$s[1]+1;
	    my $end=$s[2]+0;

	    if(exists $cov->{$chr}){
		push(@{$cov->{$chr}{"start"}},$start);
		push(@{$cov->{$chr}{"end"}},$end);
		push(@{$cov->{$chr}{"val"}},$val);
	    }
	    else{
		print "reading chr ".$chr."\n";
		$cov->{$chr}={"start"=>[$start],"end"=>[$end],"val"=>[$val]};
	    }
	}
	
	$line=<$input>;	
    }

    close($input);
}

##############################################################

sub transformCoverage{
    my $cov=$_[0];
    my $chr=$_[1];
    my $hashcov=$_[2];
        
    if(exists $cov->{$chr}){
	my $nb=@{$cov->{$chr}{"start"}};
	
	for(my $i=0; $i<$nb; $i++){
	    my $start=${$cov->{$chr}{"start"}}[$i];
	    my $end=${$cov->{$chr}{"end"}}[$i];
	    my $val=${$cov->{$chr}{"val"}}[$i];
	    
	    for(my $pos=$start; $pos<=$end; $pos++){
		$hashcov->{$pos}=$val;
	    }
	}
    }
}

##############################################################

sub transformGene{
    my $exons=$_[0];
    my $gene=$_[1];
    my $pos=$_[2];

    if(exists $exons->{$gene}){
	my $nbex=@{$exons->{$gene}{"begin"}};
	my $strand=$exons->{$gene}{"strand"};

	if($strand eq "1"){
	    for(my $i=0; $i<$nbex; $i++){
		my $begin=${$exons->{$gene}{"begin"}}[$i];
		my $end=${$exons->{$gene}{"end"}}[$i];

		for(my $k=$begin; $k<=$end; $k++){
		    push(@{$pos}, $k);
		}
	    }
	} 

	if($strand eq "-1"){
	    for(my $i=($nbex-1); $i>=0; $i--){
		my $begin=${$exons->{$gene}{"begin"}}[$i];
		my $end=${$exons->{$gene}{"end"}}[$i];

		for(my $k=$end; $k>=$begin; $k--){
		    push(@{$pos}, $k);
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
    print "This script computes read coverage variation along a gene.\n";
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
$parameters{"pathCoverageForward"}="NA";
$parameters{"pathCoverageReverse"}="NA";
$parameters{"nbWindows"}="NA";
$parameters{"minLength"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathExonBlocks","pathCoverageForward","pathCoverageReverse","nbWindows","minLength","pathOutput");
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

print "Reading exon blocks...\n";
my %exons;
readBlocks($parameters{"pathExonBlocks"},\%exons);
my $nbg=keys %exons;
print "There are ".$nbg." genes.\n";
print "Done.\n";

my %genechr;

foreach my $gene (keys %exons){
    my $chr=$exons{$gene}{"chr"};

    if(exists $genechr{$chr}){
	push(@{$genechr{$chr}}, $gene);
    }
    else{
	$genechr{$chr}=[$gene];
    }
}

my $nbgenes=keys %exons;
print "There are ".$nbgenes." genes.\n";

print "Reading coverage...\n";
my %coveragefwd;
my %coveragerev;
readCoverage($parameters{"pathCoverageForward"},\%coveragefwd);
readCoverage($parameters{"pathCoverageReverse"},\%coveragerev);
print "Done.\n";

print "Computing coverage in windows...\n";


my $nbwindows=$parameters{"nbWindows"}+0;
my $minlength=$parameters{"minLength"}+0;

print "Number of windows: ".$nbwindows."\n";
print "Minimum length: ".$minlength."\n";

open(my $output, ">".$parameters{"pathOutput"});

my $line="GeneID\tChr\tStart\tEnd\tStrand";

for(my $i=0; $i<$nbwindows; $i++){
    $line.="\tWindow".$i;
}

print $output $line."\n";


foreach my $chr (keys %genechr){
    print $chr."\n";

    my %thiscovfwd;
    if(exists $coveragefwd{$chr}){
	transformCoverage(\%coveragefwd,$chr,\%thiscovfwd);
    }
    
    my %thiscovrev;
    if(exists $coveragerev{$chr}){
	transformCoverage(\%coveragerev,$chr,\%thiscovrev);
    }

    foreach my $gene (@{$genechr{$chr}}){
	my @pos;
	transformGene(\%exons, $gene, \@pos);
	my $nbpos=@pos;
	
	my $strand=$exons{$gene}{"strand"};
	my $genestart=min @{$exons{$gene}{"begin"}};
	my $geneend=max @{$exons{$gene}{"end"}};
	
	if($nbpos>=$minlength){
	    
	    my $line=$gene."\t".$chr."\t".$genestart."\t".$geneend."\t".$strand;
	    
	    my $size=int($nbpos/$nbwindows);
	    my $reste=$nbpos-($size*$nbwindows);
	    
	    my @sizes=($size) x $nbwindows;
	    
	    for(my $i=0; $i<$reste; $i++){
		$sizes[$i]++;
	    }
	    
	    my $currentstart=0;
	    my $currentend=0;
	    
	    for(my $i=0; $i<$nbwindows; $i++){
		$currentend=$currentstart+$sizes[$i]-1;
		
		my $sumcov=0;
		
		for(my $k=$currentstart; $k<=$currentend; $k++){
		    my $thispos=$pos[$k];
		    
		    if($strand eq "1"){
			if(exists $thiscovfwd{$thispos}){
			    $sumcov+=$thiscovfwd{$thispos};
			}
		    }
			
		    if($strand eq "-1"){
			if(exists $thiscovrev{$thispos}){
			    $sumcov+=$thiscovrev{$thispos};
			}
		    }
		}
		
		my $meancov=($sumcov+0.0)/($sizes[$i]+0.0);
		
		$line.="\t".$meancov;
		
		$currentstart=$currentend+1;
	    }
	    
	    print $output $line."\n";
	}
    }
    
}

close($output);

print "Done.\n";

##############################################################
