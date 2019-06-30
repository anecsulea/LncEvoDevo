use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $refblocks=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $gene=$s[0];
	my $idexon=$s[1];
	my $chr=$s[2];
	my $begin=$s[3]+0;
	my $end=$s[4]+0;
	my $strand=$s[5];

	if(exists $refblocks->{$gene}){
	    push(@{$refblocks->{$gene}{"begin"}},$begin);
	    push(@{$refblocks->{$gene}{"end"}},$end);
	}
	else{
	    $refblocks->{$gene}={"gene"=>$gene,"chr"=>$chr,"strand"=>$strand,"begin"=>[$begin],"end"=>[$end]};
	}
	    
	$line=<$input>;
    }
    close($input);
}

##############################################################

sub extractIntrons{
    my $exons=$_[0];
    my $minsize=$_[1];
    my $maxsize=$_[2];
    my $introns=$_[3];

    foreach my $gene (keys %{$exons}){
	my $nbexons=@{$exons->{$gene}{"begin"}};
       
	if($nbexons>=2){
	    
	    my $chr=$exons->{$gene}{"chr"};
	    my $strand=$exons->{$gene}{"strand"};
	    
	    for(my $i=0; $i<($nbexons-1); $i++){
		
		my $beginintron=${$exons->{$gene}{"end"}}[$i]+1;
		my $endintron=${$exons->{$gene}{"begin"}}[$i+1]-1;

		my $size=$endintron-$beginintron+1;
		
		if($size>=$minsize && $size<=$maxsize){
		    
		    if(exists $introns->{$gene}){
			push(@{$introns->{$gene}{"begin"}},$beginintron);
			push(@{$introns->{$gene}{"end"}},$endintron);
		    }
		    else{
			$introns->{$gene}={"chr"=>$chr,"strand"=>$strand,"begin"=>[$beginintron],"end"=>[$endintron]};
		    }
		}
	    }
	}
    }
}

##############################################################

sub extractGeneCoords{
    my $refexons=$_[0];
    my $refgenes=$_[1];
    
    foreach my $gene (keys %{$refexons}){
	my $chr=$refexons->{$gene}{"chr"};
	my $strand=$refexons->{$gene}{"strand"};
	my $begin=min @{$refexons->{$gene}{"begin"}};
	my $end=max @{$refexons->{$gene}{"end"}};
	
	$refgenes->{$gene}={"chr"=>$chr,"strand"=>$strand,"begin"=>$begin,"end"=>$end};
    }
}

##############################################################

sub orderGenes{

    my $refunordered=$_[0];
    my $refordered=$_[1];
         
    my %hashbegin;

    foreach my $gene (keys %{$refunordered}){

	my $chr=$refunordered->{$gene}{"chr"};
	my $b=$refunordered->{$gene}{"begin"};
	my $e=$refunordered->{$gene}{"end"};
	
	if(exists $hashbegin{$chr}){
	    if(exists $hashbegin{$chr}{$b}){
		push(@{$hashbegin{$chr}{$b}{"end"}},$e);
		push(@{$hashbegin{$chr}{$b}{"gene"}},$gene);
		
	    }
	    else{
		$hashbegin{$chr}{$b}={"end"=>[$e],"gene"=>[$gene]};
	    }
	}
	else{
	    $hashbegin{$chr}={$b=>{"end"=>[$e],"gene"=>[$gene]}};
	}
    }
    
    foreach my $chr (keys %hashbegin){
	$refordered->{$chr}={"begin"=>[],"end"=>[],"gene"=>[]};

	my @uniquebegin=keys %{$hashbegin{$chr}}; 
	my @sortedbegin=sort {$a <=> $b} @uniquebegin;
	
	
	foreach my $begin (@sortedbegin){
	    my $nb=@{$hashbegin{$chr}{$begin}{"end"}};
	    
	    for(my $i=0;$i<$nb;$i++){
		push(@{$refordered->{$chr}{"begin"}},$begin);
		push(@{$refordered->{$chr}{"end"}},${$hashbegin{$chr}{$begin}{"end"}}[$i]);
		push(@{$refordered->{$chr}{"gene"}},${$hashbegin{$chr}{$begin}{"gene"}}[$i]);
	    }
	}
    }

}

##############################################################

sub orderIntrons{
    my $refunordered=$_[0];
    my $refordered=$_[1];
         
    my %hashbegin;

    foreach my $gene (keys %{$refunordered}){
	my $chr=$refunordered->{$gene}{"chr"};
	my $strand=$refunordered->{$gene}{"strand"};
	my $nb=@{$refunordered->{$gene}{"begin"}};
	
	for(my $i=0;$i<$nb;$i++){

	    my $b=${$refunordered->{$gene}{"begin"}}[$i];
	    my $e=${$refunordered->{$gene}{"end"}}[$i];
	   

	     my $key=$gene.",".$chr.",".$b.",".$e.",".$strand;
 
	    if(exists $hashbegin{$chr}){
		if(exists $hashbegin{$chr}{$b}){
		    push(@{$hashbegin{$chr}{$b}{"end"}},$e);
		    push(@{$hashbegin{$chr}{$b}{"key"}},$key);
		}
		else{
		    $hashbegin{$chr}{$b}={"end"=>[$e],"key"=>[$key]};
		}
	    }
	    else{
		$hashbegin{$chr}={$b=>{"end"=>[$e],"key"=>[$key]}};
	    }
	}
    }
    
    foreach my $chr (keys %hashbegin){
	$refordered->{$chr}={"begin"=>[],"end"=>[],"key"=>[]};

	my @uniquebegin=keys %{$hashbegin{$chr}}; 
	my @sortedbegin=sort {$a <=> $b} @uniquebegin;
	
	
	foreach my $begin (@sortedbegin){
	    my $nb=@{$hashbegin{$chr}{$begin}{"end"}};
	    
	    for(my $i=0;$i<$nb;$i++){
		push(@{$refordered->{$chr}{"begin"}},$begin);
		push(@{$refordered->{$chr}{"end"}},${$hashbegin{$chr}{$begin}{"end"}}[$i]);
		push(@{$refordered->{$chr}{"key"}},${$hashbegin{$chr}{$begin}{"key"}}[$i]);
	    }
	}
    }

}

##############################################################

sub checkOverlappingIntrons{
    ## look for introns that overlap with other genes
    
    my $reforderedintrons=$_[0];
    my $reforderedgenes=$_[1];
    my $refoverlap=$_[2];
    
    my %removed;

    foreach my $chr (keys %{$reforderedintrons}){
	my $nbintrons=@{$reforderedintrons->{$chr}{"begin"}};
	my $nbgenes=@{$reforderedgenes->{$chr}{"begin"}};

	my $firstgene=0;
	
	for(my $i=0;$i<$nbintrons;$i++){
	    my $begin=${$reforderedintrons->{$chr}{"begin"}}[$i];
	    my $end=${$reforderedintrons->{$chr}{"end"}}[$i];
	    my $key=${$reforderedintrons->{$chr}{"key"}}[$i];
	    my @s=split(",",$key);
	    my $gene=$s[0];
	    
	    my $j=$firstgene;

	    while($j<$nbgenes && ${$reforderedgenes->{$chr}{"end"}}[$j]<$begin){
		$j++;
	    }
	    
	    $firstgene=$j;
	    
	    while($j<$nbgenes && ${$reforderedgenes->{$chr}{"begin"}}[$j]<=$end){
		my $thisgene = ${$reforderedgenes->{$chr}{"gene"}}[$j];
		my $thisbegin = ${$reforderedgenes->{$chr}{"begin"}}[$j];
		my $thisend = ${$reforderedgenes->{$chr}{"end"}}[$j];
	       
		my $M=max($thisbegin,$begin);
		my $m=min($thisend,$end);

		if($M<=$m){
		    if(!($thisgene eq $gene)){
			if(exists $refoverlap->{$key}){
			    $refoverlap->{$key}{$thisgene}=1;
			}
			else{
			    $refoverlap->{$key}={$thisgene=>1};
			}
		    }
		}
		
		$j++;
	    }
	}

    }

}

##############################################################
##############################################################

my %parameters;
$parameters{"pathExonBlocks"}="NA";
$parameters{"minIntronSize"}="NA";
$parameters{"maxIntronSize"}="NA";
$parameters{"excludeLength"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathExonBlocks","minIntronSize","maxIntronSize","excludeLength","pathOutput");

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

print "Reading exon blocks...\n";
my %exons;
readExonBlocks($parameters{"pathExonBlocks"},\%exons);
my $nbgenes=keys %exons;
print "There are ".$nbgenes." genes.\n";
print "Done.\n";

print "Extracting and ordering genes...\n";
my %genes;
extractGeneCoords(\%exons,\%genes);
my %orderedgenes;
orderGenes(\%genes,\%orderedgenes);
print "Done.\n";

print "Extracting and ordering introns...\n";
my %introns;
my $minsize=$parameters{"minIntronSize"}+0;
my $maxsize=$parameters{"maxIntronSize"}+0;
my $exLen=$parameters{"excludeLength"}+0;

extractIntrons(\%exons,$minsize,$maxsize,\%introns);
my %orderedintrons;
orderIntrons(\%introns, \%orderedintrons);
print "Done.\n";

print "Checking overlap between genes and introns...\n";
my %overlap;
checkOverlappingIntrons(\%orderedintrons,\%orderedgenes,\%overlap);
my $nbov=keys %overlap;
print "There are ".$nbov." introns that overlap with other genes.\n";
print "Done.\n";

print "Writing output...\n";
open(my $output,">".$parameters{"pathOutput"});
foreach my $gene (keys %introns){
    my $chr=$introns{$gene}{"chr"};
    my $strand=$introns{$gene}{"strand"};
    
    my $nbi=@{$introns{$gene}{"begin"}};

    for(my $i=0; $i<$nbi; $i++){
	my $begin=${$introns{$gene}{"begin"}}[$i];
	my $end=${$introns{$gene}{"end"}}[$i];
	my $key=$gene.",".$chr.",".$begin.",".$end.",".$strand;
	
	if(!(exists $overlap{$key})){
	    my $startregion=$begin+$exLen;
	    my $endregion=$end-$exLen;
	    print $output $gene."\t".$gene.".".$i."\t".$chr."\t".$startregion."\t".$endregion."\t".$strand."\t".$begin."\t".$end."\n";
	}
    }
}
close($output);
print "Done.\n";

##############################################################
