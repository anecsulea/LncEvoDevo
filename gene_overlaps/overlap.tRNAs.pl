use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
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
	my $start=$s[3]+0;
	my $end=$s[4]+0;
	my $strand=$s[5];

	if($strand eq "."){
	    $strand="NA";
	}

	if(exists $refblocks->{$gene}){
	    push(@{$refblocks->{$gene}{"start"}},$start);
	    push(@{$refblocks->{$gene}{"end"}},$end);
	}
	else{
	    $refblocks->{$gene}={"chr"=>$chr,"strand"=>$strand,"start"=>[$start],"end"=>[$end]};
	}
	    
	$line=<$input>;
    }
    close($input);
}

##########################################################################

sub orderExonBlocks{
    my $refblocks=$_[0];
    my $refordered=$_[1];
    
    my %hashstart;
    
    foreach my $gene (keys %{$refblocks}){
	my $nbblocks=@{$refblocks->{$gene}{"start"}};
	my $chr=$refblocks->{$gene}{"chr"};
	my $strand=$refblocks->{$gene}{"strand"};

	for(my $i=0;$i<$nbblocks;$i++){
	    my $tb=${$refblocks->{$gene}{"start"}}[$i];
	    my $te=${$refblocks->{$gene}{"end"}}[$i];
	
	    if(exists $hashstart{$chr}){
		if(exists $hashstart{$chr}{$tb}){
		    if(exists $hashstart{$chr}{$tb}{$te}){
			if(exists $hashstart{$chr}{$tb}{$te}{$strand}){
			    push(@{$hashstart{$chr}{$tb}{$te}{$strand}},$gene);
			}
			else{
			    $hashstart{$chr}{$tb}{$te}{$strand}=[$gene];
			}
		    }
		    else{
			$hashstart{$chr}{$tb}{$te}={$strand=>[$gene]};
		    }
		}
		else{
		    $hashstart{$chr}{$tb}={$te=>{$strand=>[$gene]}};
		}
	    }
	    else{
		$hashstart{$chr}={$tb=>{$te=>{$strand=>[$gene]}}};
	    }
	}
	
    }

    foreach my $chr (keys %hashstart){
	$refordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"gene"=>[]};

	my @uniquestart=keys %{$hashstart{$chr}};
	
	my @sortedstart = sort {$a <=> $b} @uniquestart;

	foreach my $start (@sortedstart){
	    my @uniqueend=keys %{$hashstart{$chr}{$start}};
	    my @sortedend = sort {$a <=> $b} @uniqueend;

	    foreach my $end (@sortedend){
		foreach my $strand (keys %{$hashstart{$chr}{$start}{$end}}){
		    foreach my $gene (@{$hashstart{$chr}{$start}{$end}{$strand}}){
			push(@{$refordered->{$chr}{"start"}},$start);
			push(@{$refordered->{$chr}{"end"}},$end);
			push(@{$refordered->{$chr}{"strand"}},$strand);
			push(@{$refordered->{$chr}{"gene"}},$gene);
		    }
		} 
	    }
	} 
    }

}

##############################################################

sub orderCoords{
    my $exons=$_[0];
    my $refordered=$_[1];

    my %hashstart;
    
    foreach my $exid (keys %{$exons}){
	my $chr=$exons->{$exid}{"chr"};
	my $b=$exons->{$exid}{"start"};
	my $e=$exons->{$exid}{"end"};
	my $s=$exons->{$exid}{"strand"};
	
	if(exists $hashstart{$chr}){
	    if(exists $hashstart{$chr}{$b}){
		push(@{$hashstart{$chr}{$b}{"end"}},$e);
		push(@{$hashstart{$chr}{$b}{"strand"}},$s);
		push(@{$hashstart{$chr}{$b}{"id"}},$exid);
	    }
	    else{
		$hashstart{$chr}{$b}={"end"=>[$e],"strand"=>[$s],"id"=>[$exid]};
	    }
	}
	else{
	    $hashstart{$chr}={$b=>{"end"=>[$e],"strand"=>[$s],"id"=>[$exid]}};
	}
    }
    
    foreach my $chr (keys %hashstart){
	$refordered->{$chr}={"start"=>[], "end"=>[], "id"=>[], "strand"=>[]};
	
	my @uniquestart=keys %{$hashstart{$chr}};
	
	my @sortedstart = sort {$a <=> $b} @uniquestart;
	
	foreach my $b (@sortedstart){
	    
	    my $nbblocks=@{$hashstart{$chr}{$b}{"end"}};
	    
	    for(my $i=0;$i<$nbblocks;$i++){
		my $strand=${$hashstart{$chr}{$b}{"strand"}}[$i];
		
		push(@{$refordered->{$chr}{"start"}},$b);
		push(@{$refordered->{$chr}{"strand"}},$strand);
		push(@{$refordered->{$chr}{"end"}},${$hashstart{$chr}{$b}{"end"}}[$i]);
		push(@{$refordered->{$chr}{"id"}},${$hashstart{$chr}{$b}{"id"}}[$i]);
	    }
	}
    }
}

##############################################################

sub readTRNAs{
    my $pathin=$_[0];
    my $pirnas=$_[1];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    my $line=<$input>; ## header
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $chr=$s[1];
	my $prefix=substr $chr, 0,3;
	
	if($prefix eq "chr"){
	    $chr=substr $chr,3;
	}
	
	my $start=$s[2]+1; ## 1-based
	my $end=$s[3]+0;
	
	my $ss=$s[6];
	my $strand="NA";
	
	if($ss eq "+"){
	    $strand="1";
	} else{
	    if($ss eq "-"){
		$strand="-1";
	    } else{
		print "Weird strand!\n";
		print $line."\n";
		exit(1);
	    } 
	}
	
	my $id=$chr.",".$start.",".$end.",".$strand;
	
	$pirnas->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	
	$line=<$input>;
    }
    
    close($input);

}

##############################################################

sub overlapBlocks{
    my $refblocks1=$_[0]; ## ordered exon blocks
    my $refblocks2=$_[1]; ## ordered exon blocks
    my $type=$_[2];
    my $refoverlap=$_[3]; ## overlap with blocks
   
    foreach my $chr (keys %{$refblocks1}){
	if(exists $refblocks2->{$chr}){
	    my $nbblocks1=@{$refblocks1->{$chr}{"start"}};
	    my $nbblocks2=@{$refblocks2->{$chr}{"start"}};
	    
	    my $firstindex2=0;  ## this is where we start looking for overlap
	    
	    for(my $i=0; $i<$nbblocks1; $i++){
		my $start1=${$refblocks1->{$chr}{"start"}}[$i];
		my $end1=${$refblocks1->{$chr}{"end"}}[$i];
		my $gene1=${$refblocks1->{$chr}{"gene"}}[$i];
		my $strand1=${$refblocks1->{$chr}{"strand"}}[$i];

		my $j=$firstindex2;

		while($j<$nbblocks2 && ${$refblocks2->{$chr}{"end"}}[$j]<$start1){ ## there cannnot be any overlap before that 
		    $j++;
		}

		$firstindex2=$j;

		while($j<$nbblocks2 && ${$refblocks2->{$chr}{"start"}}[$j]<=$end1){  ## we stop looking for overlap if the start coordinate of the second set of blocks is larger than the end coordinate of block
		 
		    my $M=max($start1,${$refblocks2->{$chr}{"start"}}[$j]);
		    my $m=min($end1,${$refblocks2->{$chr}{"end"}}[$j]);
		    
		    my $strand2=${$refblocks2->{$chr}{"strand"}}[$j];

		    if(($type eq "sense" && $strand1 eq $strand2) || ($type eq "antisense" && $strand1 ne $strand2) || ($type eq "any")){
			my $lenov=$m-$M+1;
			
			if($lenov>=1){
			    if(exists $refoverlap->{$gene1}){
				push(@{$refoverlap->{$gene1}{"start"}},$M);
				push(@{$refoverlap->{$gene1}{"end"}},$m);
			    }
			    else{
				$refoverlap->{$gene1}={"start"=>[$M],"end"=>[$m]};
			    }
			}
		    }
		    
		    $j++;
		}
	    }
	}
    }
}

##############################################################

sub makeBlocks{
    my $intervals=$_[0];
    my $blocks=$_[1];

    my %hashpos;

    my $nb=@{$intervals->{"start"}};
    $blocks->{"start"}=[];
    $blocks->{"end"}=[];

    for(my $i=0; $i<$nb; $i++){
	my $start=${$intervals->{"start"}}[$i];
	my $end=${$intervals->{"end"}}[$i];

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
    my $totlen=0;

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
	    push(@{$blocks->{"start"}}, $currentstart);
	    push(@{$blocks->{"end"}}, $currentend);
	    $totlen+=($currentend-$currentstart+1);

	    $currentstart=$thisstart;
	    $currentend=$thisend;
	}
    }
    
    push(@{$blocks->{"start"}}, $currentstart);
    push(@{$blocks->{"end"}}, $currentend);
    $totlen+=($currentend-$currentstart+1);

    $blocks->{"totallength"}=$totlen;
    
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes overlap with tRNAs.\n";
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
$parameters{"pathTRNAs"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathExonBlocks", "pathTRNAs", "pathOutput");

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

print "Reading annotations...\n";

my %exonblocks;

readExonBlocks($parameters{"pathExonBlocks"}, \%exonblocks);

my $nbg=keys %exonblocks;

print "Found ".$nbg." genes.\n";

print "Done.\n";

##############################################################

print "Reading tRNAs...\n";

my %trnas;
readTRNAs($parameters{"pathTRNAs"},  \%trnas);

my $nbpi=keys %trnas;

print "Found ".$nbpi." tRNAs.\n";

print "Done.\n";

##############################################################

print "Ordering exons...\n";

my %orderedexons;

orderExonBlocks(\%exonblocks, \%orderedexons);

my %orderedtrnas;

orderCoords(\%trnas, \%orderedtrnas);

print "Done.\n";

##############################################################

print "Computing overlap with tRNAs...\n";

my %overlapsense;

overlapBlocks(\%orderedexons, \%orderedtrnas, "sense", \%overlapsense);

my $nbov=keys %overlapsense;

print "Found ".$nbov." genes with tRNA coverage on the sense strand.\n";


my %overlapantisense;

overlapBlocks(\%orderedexons, \%orderedtrnas, "antisense", \%overlapantisense);

my $nbov=keys %overlapantisense;

print "Found ".$nbov." genes with tRNA coverage on the antisense strand.\n";


my %overlapboth;

overlapBlocks(\%orderedexons, \%orderedtrnas, "any", \%overlapboth);

my $nbov=keys %overlapboth;

print "Found ".$nbov." genes with tRNA coverage on any strand.\n";
   
print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tTotalExonicLength\tLengthOverlapSense\tLengthOverlapAntisense\tLengthOverlapBothStrands\n";

foreach my $gene (keys %exonblocks){
    my $totlen=0;
    
    my $nbex=@{$exonblocks{$gene}{"start"}};
    
    for(my $i=0; $i<$nbex; $i++){
	$totlen+=(${$exonblocks{$gene}{"end"}}[$i] - ${$exonblocks{$gene}{"start"}}[$i]+1);
    }

    my $lensense=0;
    my $lenantisense=0;
    my $lenboth=0;

    if(exists $overlapsense{$gene}){
	my %blocks;
	makeBlocks($overlapsense{$gene}, \%blocks);
	my $nbb=@{$blocks{"start"}};

	for(my $i=0; $i<$nbb; $i++){
	    $lensense+=(${$blocks{"end"}}[$i] - ${$blocks{"start"}}[$i]+1);
	}
    }

    if(exists $overlapantisense{$gene}){
	my %blocks;
	makeBlocks($overlapantisense{$gene}, \%blocks);
	my $nbb=@{$blocks{"start"}};

	for(my $i=0; $i<$nbb; $i++){
	    $lenantisense+=(${$blocks{"end"}}[$i] - ${$blocks{"start"}}[$i]+1);
	}
    }

    if(exists $overlapboth{$gene}){
	my %blocks;
	makeBlocks($overlapboth{$gene}, \%blocks);
	my $nbb=@{$blocks{"start"}};
	
	for(my $i=0; $i<$nbb; $i++){
	    $lenboth+=(${$blocks{"end"}}[$i] - ${$blocks{"start"}}[$i]+1);
	}
    }

    print $output $gene."\t".$totlen."\t".$lensense."\t".$lenantisense."\t".$lenboth."\n";
}

close($output);

print "Done.\n";

##############################################################
