use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readNonoverlappingAnnot{
    my $pathin=$_[0];
    my $exons=$_[1];
    my $genestx=$_[2];
    my $txexons=$_[3];

    open(my $input, $pathin);
    
    my $line=<$input>;
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
	
	my $start=$s[$header{"NonoverlapStart"}]+0;
	my $end=$s[$header{"NonoverlapEnd"}]+0;
	
	if($start>0 && $end>0){
	    my $geneid=$s[$header{"GeneID"}];
	    my $txid=$s[$header{"TranscriptID"}];
	    my $chr=$s[$header{"Chr"}];
	    my $strand=$s[$header{"Strand"}];
	    
	    my $exonid=$chr.",".$start.",".$end.",".$strand;
	  
	    if(exists $exons->{$exonid}){
		$exons->{$exonid}{"transcripts"}{$txid}=1;
		$exons->{$exonid}{"genes"}{$geneid}=1;
	    } else{
		$exons->{$exonid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand, "genes"=>{$geneid=>1}, "transcripts"=>{$txid=>1}};
	    }

	    if(exists $genestx->{$geneid}){
		$genestx->{$geneid}{$txid}=1;
	    }
	    else{
		$genestx->{$geneid}={$txid=>1};
	    }

	    if(exists $txexons->{$txid}){
		$txexons->{$txid}{$exonid}=1;
	    }
	    else{
		$txexons->{$txid}={$exonid=>1};
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

sub orderExons{
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
	$refordered->{$chr}={};
	
	my @uniquestart=keys %{$hashstart{$chr}};
	
	my @sortedstart = sort {$a <=> $b} @uniquestart;
	
	foreach my $b (@sortedstart){
	    
	    my $nbblocks=@{$hashstart{$chr}{$b}{"end"}};
	    
	    for(my $i=0;$i<$nbblocks;$i++){
		my $strand=${$hashstart{$chr}{$b}{"strand"}}[$i];
		
		if(!(exists $refordered->{$chr}{$strand})){
		    $refordered->{$chr}{$strand}={"start"=>[], "end"=>[], "id"=>[]};
		}
		
		push(@{$refordered->{$chr}{$strand}{"start"}},$b);
		push(@{$refordered->{$chr}{$strand}{"end"}},${$hashstart{$chr}{$b}{"end"}}[$i]);
		push(@{$refordered->{$chr}{$strand}{"id"}},${$hashstart{$chr}{$b}{"id"}}[$i]);
	    }
	}
    }
}

##############################################################

sub readCSFCoveredRegions{
    my $pathin=$_[0];
    my $regions=$_[1];

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
	my @s=split("\t" , $line);
	
	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	
	my $prefix=substr $chr,0,3;
	
	if($prefix eq "chr"){
	    $chr=substr $chr,3;
	}

	if($strand eq "+"){
	    $strand="1";
	}
	
	if($strand eq "-"){
	    $strand="-1";
	}
	
	if(exists $regions->{$chr}){
	    if(exists $regions->{$chr}{$strand}){
		my $laststart=${$regions->{$chr}{$strand}{"start"}}[-1];

		## we can have windows at exactly the same position if the multiple genome alignment was done on a different reference, we'll use the maximum score

		if($start<$laststart){
		    print "Data are not ordered! ".$start." ".$laststart."\n";
		    exit(1);
		}

		push(@{$regions->{$chr}{$strand}{"start"}}, $start);
		push(@{$regions->{$chr}{$strand}{"end"}}, $end);
	    }
	    else{
		$regions->{$chr}{$strand}={"start"=>[$start], "end"=>[$end]};
	    }
	}
	else{
	    $regions->{$chr}={$strand=>{"start"=>[$start], "end"=>[$end]}};
	}
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub readPositiveCSFScores{
    my $pathin=$_[0];
    my $regions=$_[1];

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
	my @s=split("\t" , $line);
	
	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	my $sumscore=$s[$header{"SumScore"}]+0.0;
	my $nbinfo=$s[$header{"NbInfoCodons"}]+0;

	my $prefix=substr $chr,0,3;

	if($prefix eq "chr"){
	    $chr=substr $chr,3;
	}
	
	if($strand eq "+"){
	    $strand="1";
	}

	if($strand eq "-"){
	    $strand="-1";
	}

	my $meanscore=$sumscore/$nbinfo;

	if(exists $regions->{$chr}){
	    if(exists $regions->{$chr}{$strand}){
		my $laststart=${$regions->{$chr}{$strand}{"start"}}[-1];
		
		## we can have windows at the exactly same position if the multiple genome alignment was done on a different reference, we'll use the maximum score
		
		if($start<$laststart){ 
		    print "Data are not ordered! ".$start." ".$laststart."\n";
		    exit(1);
		}

		push(@{$regions->{$chr}{$strand}{"start"}}, $start);
		push(@{$regions->{$chr}{$strand}{"end"}}, $end);
		push(@{$regions->{$chr}{$strand}{"meanscore"}}, $meanscore);
	    }
	    else{
		$regions->{$chr}{$strand}={"start"=>[$start], "end"=>[$end], "meanscore"=>[$meanscore]};
	    }
	}
	else{
	    $regions->{$chr}={$strand=>{"start"=>[$start], "end"=>[$end], "meanscore"=>[$meanscore]}};
	}
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub computeCSFCoverage{
    my $orderedexons=$_[0];
    my $regions=$_[1];
    my $coverage=$_[2];

    foreach my $chr (keys %{$orderedexons}){
	foreach my $strand (keys %{$orderedexons->{$chr}}){
	   	      
	    if(exists $regions->{$chr}){
		if(exists $regions->{$chr}{$strand}){
		    my $nbex=@{$orderedexons->{$chr}{$strand}{"start"}};
		    my $nbreg=@{$regions->{$chr}{$strand}{"start"}};
		    
		    my $firstreg=0;

		    for(my $i=0; $i<$nbex; $i++){
			my $startex=${$orderedexons->{$chr}{$strand}{"start"}}[$i];
			my $endex=${$orderedexons->{$chr}{$strand}{"end"}}[$i];
			my $idex=${$orderedexons->{$chr}{$strand}{"id"}}[$i];

			my $j=$firstreg;

			while($j<$nbreg && ${$regions->{$chr}{$strand}{"end"}}[$j] < $startex){
			    $j++;
			}

			$firstreg=$j;

			while($j<$nbreg && ${$regions->{$chr}{$strand}{"start"}}[$j] <= $endex){
			    my $startreg=${$regions->{$chr}{$strand}{"start"}}[$j];
			    my $endreg=${$regions->{$chr}{$strand}{"end"}}[$j];
			    
			    my $M=max($startex, $startreg);
			    my $m=min($endex, $endreg);

			    if($M<=$m){
				if(exists $coverage->{$idex}){
				    push(@{$coverage->{$idex}{"start"}}, $M);
				    push(@{$coverage->{$idex}{"end"}}, $m);
				}
				else{
				    $coverage->{$idex}={"start"=>[$M], "end"=>[$m]};
				}
			    }
			    
			    $j++;
			}
		    }
		}
	    }
	}
    }	
}

##############################################################

sub overlapPositiveScores{
    my $orderedexons=$_[0];
    my $positivescores=$_[1];
    my $minoverlap=$_[2];
    my $overlap=$_[3];

    foreach my $chr (keys %{$orderedexons}){
	foreach my $strand (keys %{$orderedexons->{$chr}}){

	    if(exists $positivescores->{$chr}){
		if(exists $positivescores->{$chr}{$strand}){
		    my $nbex=@{$orderedexons->{$chr}{$strand}{"start"}};
		    my $nbreg=@{$positivescores->{$chr}{$strand}{"start"}};
		    
		    my $firstreg=0;

		    for(my $i=0; $i<$nbex; $i++){
			my $startex=${$orderedexons->{$chr}{$strand}{"start"}}[$i];
			my $endex=${$orderedexons->{$chr}{$strand}{"end"}}[$i];
			my $idex=${$orderedexons->{$chr}{$strand}{"id"}}[$i];

			my $j=$firstreg;

			while($j<$nbreg && ${$positivescores->{$chr}{$strand}{"end"}}[$j] < $startex){
			    $j++;
			}

			$firstreg=$j;

			while($j<$nbreg && ${$positivescores->{$chr}{$strand}{"start"}}[$j] <= $endex){
			    my $startreg=${$positivescores->{$chr}{$strand}{"start"}}[$j];
			    my $endreg=${$positivescores->{$chr}{$strand}{"end"}}[$j];
			    
			    my $M=max($startex, $startreg);
			    my $m=min($endex, $endreg);

			    my $lenov=($m-$M+1);

			    if($lenov>=$minoverlap){
				my $score=${$positivescores->{$chr}{$strand}{"meanscore"}}[$j];

				if(exists $overlap->{$idex}){
				    push(@{$overlap->{$idex}{"start"}}, $M);
				    push(@{$overlap->{$idex}{"end"}}, $m);
				    push(@{$overlap->{$idex}{"score"}}, $score);
				}
				else{
				    $overlap->{$idex}={"start"=>[$M], "end"=>[$m], "score"=>[$score]};
				}
			    }
			    
			    $j++;
			}
		    }
		}
	    }
	}	
    }
}

##########################################################################

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

    my $nbsorted=@sortedpos;

    if($nbsorted>0){
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
    }

    $blocks->{"totallength"}=$totlen;
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes CSF statistics for exons and genes.\n";
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
$parameters{"pathNonOverlappingAnnot"}="NA";
$parameters{"pathCoveredRegions"}="NA";
$parameters{"pathPositiveScores"}="NA";
$parameters{"minFractionOverlap"}="NA";
$parameters{"minLengthOverlap"}="NA";
$parameters{"pathOutputExons"}="NA";
$parameters{"pathOutputTranscripts"}="NA";
$parameters{"pathOutputGenes"}="NA";

my @defaultpars=("pathNonOverlappingAnnot", "pathCoveredRegions", "pathPositiveScores", "minFractionOverlap", "minLengthOverlap", "pathOutputExons", "pathOutputTranscripts", "pathOutputGenes");

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

my %exons;
my %genestx;
my %txexons;

readNonoverlappingAnnot($parameters{"pathNonOverlappingAnnot"}, \%exons, \%genestx, \%txexons);

my $nbex=keys %exons;
my $nbg=keys %genestx;
my $nbtx=keys %txexons;

print "Found ".$nbex." exons, ".$nbtx." transcripts and ".$nbg." genes.\n";

print "Done.\n";

##############################################################


print "Reading regions covered by CSF scores...\n";

my %coveredregions;

readCSFCoveredRegions($parameters{"pathCoveredRegions"}, \%coveredregions);

print "Done.\n";

##############################################################

print "Reading positive CSF scores...\n";

my %positivescores;

readPositiveCSFScores($parameters{"pathPositiveScores"}, \%positivescores);

print "Done.\n";

##############################################################

print "Ordering exons...\n";

my %orderedexons;

orderExons(\%exons, \%orderedexons);

print "Done.\n";

##############################################################

print "Computing overlap with covered regions...\n";

my %exoncoverage;

computeCSFCoverage(\%orderedexons, \%coveredregions, \%exoncoverage);

my $nbov=keys %exoncoverage;

print "Found ".$nbov." exons with CSF coverage.\n";
   
print "Done.\n";

##############################################################

print "Computing overlap with positive scores...\n";

my %exonscores;

overlapPositiveScores(\%orderedexons, \%positivescores, 1, \%exonscores);

my $nbov=keys %exonscores;

print "Found ".$nbov." exons with positive scores.\n";
 
print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutputExons"});

print $output "ExonID\tChr\tStart\tEnd\tStrand\tGeneID\tTranscriptID\tTotalLength\tLengthPositiveCSF\tMaxScore\tLengthCoveredCSF\n";

my %blockscoverage;
my %blocksscores;

foreach my $exonid (keys %exons){
    my $geneid=join(";", keys %{$exons{$exonid}{"genes"}});
    my $txid=join(";", keys %{$exons{$exonid}{"transcripts"}});
    my $chr=$exons{$exonid}{"chr"};
    my $start=$exons{$exonid}{"start"};
    my $end=$exons{$exonid}{"end"};
    my $strand=$exons{$exonid}{"strand"};

    my $totlen=$end-$start+1;

    my $lenpos=0;
    my $lencov=0;

    my $maxscore=0;

    if(exists $exonscores{$exonid}){
	$blocksscores{$exonid}={};

	$maxscore=max(@{$exonscores{$exonid}{"score"}});

	makeBlocks($exonscores{$exonid}, $blocksscores{$exonid});
	$lenpos=$blocksscores{$exonid}{"totallength"};
    }
    
    if(exists $exoncoverage{$exonid}){
	$blockscoverage{$exonid}={};
	
	makeBlocks($exoncoverage{$exonid}, $blockscoverage{$exonid});
	$lencov=$blockscoverage{$exonid}{"totallength"};
    }
    
    print $output $exonid."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$geneid."\t".$txid."\t".$totlen."\t".$lenpos."\t".$maxscore."\t".$lencov."\n";
}

close($output);

print "Done.\n";

##############################################################

print "Writing output and classification for genes...\n";

open(my $outputgenes, ">".$parameters{"pathOutputGenes"});
open(my $outputtx, ">".$parameters{"pathOutputTranscripts"});

my $minfroverlap=$parameters{"minFractionOverlap"}+0.0;
my $minlenoverlap=$parameters{"minLengthOverlap"}+0;

print "Minimum length overlap: ".$minlenoverlap." minimum fraction overlap: ".$minfroverlap."\n";

print $outputgenes "GeneID\tClass\tCodingTranscripts\n";
print $outputtx "GeneID\tTranscriptID\tTotalLength\tLengthCoveredCSF\tLengthOverlapPositiveCSF\tTranscriptClass\n";

foreach my $gene (keys %genestx){
    my $classgene="NA";
    my %codingtranscripts;
    my %noncodingtranscripts;
 
    foreach my $tx (keys %{$genestx{$gene}}){
	my $totlen=0;
	my $lenpositive=0;
	my $lencovered=0;
	
	foreach my $exonid (keys %{$txexons{$tx}}){
	    my $start=$exons{$exonid}{"start"};
	    my $end=$exons{$exonid}{"end"};
	    
	    $totlen+=($end-$start+1);
	    
	    if(exists $blockscoverage{$exonid}){
		$lencovered+=$blockscoverage{$exonid}{"totallength"};
	    }

	    if(exists $blocksscores{$exonid}){
		$lenpositive+=$blocksscores{$exonid}{"totallength"};
	    }
	}


	my $classtx="NA";

	my $frpos=0;
	
	if($lencovered>0){ ## the transcript has to be covered by CSF scores
	    $frpos=($lenpositive+0.0)/($totlen+0.0);
	    
	    if($frpos>=$minfroverlap  && $lenpositive>=$minlenoverlap){
		$classtx="coding";
		$classgene="coding";
		$codingtranscripts{$tx}=1;
	    }
	    else{
		$classtx="noncoding";
		$noncodingtranscripts{$tx}=1;
	    }
	}	
	
	print $outputtx $gene."\t".$tx."\t".$totlen."\t".$lencovered."\t".$lenpositive."\t".$classtx."\n";
    }

    my $idtx="NA";
    my $nbcoding=keys %codingtranscripts;
    my $nbnoncoding=keys %noncodingtranscripts;
    
    if($nbcoding>0){
	$idtx=join(",",sort(keys %codingtranscripts));
    }
    else{
	if($nbnoncoding>0){
	    $classgene="noncoding";
	}
    }
    
    print $outputgenes $gene."\t".$classgene."\t".$idtx."\n";
}

close($outputgenes);

print "Done.\n";

##############################################################
