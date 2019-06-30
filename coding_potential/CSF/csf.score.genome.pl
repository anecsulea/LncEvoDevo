use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readAlignmentsCSF{
    my $pathin=$_[0];
    my $refsp=$_[1];
    my $selectedinformants=$_[2];
    my $codingmat=$_[3];
    my $noncodingmat=$_[4];
    my $geneticcode=$_[5];
    my $penalty=$_[6];
    my $mininformants=$_[7];  
    my $startregion=$_[8];
    my $endregion=$_[9];
    my $codoncsf=$_[10];
    my $forbidden=$_[11];
    my $chrsizesref=$_[12];

    my $input;

    my @s=split("\\.",$pathin);
    my $ext=$s[-1];

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
   
    my %aln;

    while($line){

	chomp $line;
	$firstchar=substr $line,0,1;
	
	if($firstchar eq "a"){
	    
	    ## we start a new alignment

	    if($indexaln>0){
		## we already have an alignment
		## check if this alignment has the reference species

		if($indexaln%10000==0){
		    print "done ".$indexaln." alignments\n";
		}
		
		my $nbsp=keys %aln;
		
		if((exists $aln{$refsp}) && ($nbsp>=($mininformants+1))){

		    my $lenaln=length $aln{$refsp}{"sequence"};
		    
		    if($lenaln>=3){
			my $thisstart=$aln{$refsp}{"start"};
			my $thissize=length $aln{$refsp}{"sequence"};
			
			my $pos1=$thisstart-$thissize-1;
			my $pos2=$thisstart+$thissize+1;
			
			my $M=max($pos1,$startregion);
			my $m=min($pos2,$endregion);

			if($M<=$m){
			    my %nogaps;
			    removeGapsRef(\%aln,$refsp,\%nogaps);
			    computeCSFCodons(\%nogaps, $refsp, $codingmat, $noncodingmat, $geneticcode, $mininformants,$penalty, $codoncsf,$forbidden);
			}
		    }
		}
	    }

	    $indexaln++;

	    my @s=split(" ",$line);
	    my $score=$s[1];
	    my @t=split("=",$score);
	    $score=$t[1]+0;
	    $currentscore=$score;
	    %aln=();
	}

	if($firstchar eq "s"){
	    my @s=split(" ",$line);
	    my @t=split("\\.",$s[1]);
	    my $idsp=$t[0];

	    if(exists $selectedinformants->{$idsp} || $idsp eq $refsp){
		
		my $chr=$t[1];
		my $start=$s[2]+0;
		my $strand=$s[4];
		my $chrsize=$s[5]+0;	
		my $sequence=$s[6];

		if($idsp eq $refsp){
		    $chrsizesref->{$chr}=$chrsize;
		}

		if($strand eq "+"){
		    $start=$start+1; ## positions start at 1, like for our annotations
		}
		else{
		    if($strand eq "-"){
			$start=$chrsize-$start; ## position start at 1, increment will be negative
		    }
		    else{
			print "Weird strand in line ".$line."\n";
			exit(1);
		    }
		}
		
		$aln{$idsp}={"chr"=>$chr,"start"=>$start,"strand"=>$strand,"sequence"=>$sequence};
	    }
	}
		
	$line=<$input>;
    }

    close($input);
     
    ## don't forget the last alignment
    
    my $nbsp=keys %aln;
    
    if((exists $aln{$refsp}) && ($nbsp>=($mininformants+1))){
	 my $lenaln=length $aln{$refsp}{"sequence"};
	 
	 if($lenaln>=3){

	     	my $thisstart=$aln{$refsp}{"start"};
		my $thissize=length $aln{$refsp}{"sequence"};
		
		my $pos1=$thisstart-$thissize-1;
		my $pos2=$thisstart+$thissize+1;
		
		my $M=max($pos1,$startregion);
		my $m=min($pos2,$endregion);
		
		if($M<=$m){
		    my %nogaps;
		    removeGapsRef(\%aln,$refsp,\%nogaps);
		    computeCSFCodons(\%nogaps, $refsp, $codingmat, $noncodingmat, $geneticcode, $mininformants,$penalty, $codoncsf,$forbidden);
		}
	 }
    }
}

##############################################################

sub readGeneticCode{
    my $pathin=$_[0];
    my $code=$_[1];

    open(my $input,$pathin);
    my $line=<$input>;

    while($line){
	my @s=split(" ",$line);
	$code->{$s[0]}=$s[1];

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readMatrices{
    my $pathin=$_[0];
    my $matrices=$_[1];
    
    open(my $input, $pathin);

    my $sp1="NA";
    my $sp2="NA";

    my $line=<$input>;
    
    while($line){
	chomp $line;
	my $prefix=substr $line,0,1;
	
	if($prefix eq "#"){
	    $line=substr $line,2;
	    my @s=split(" ",$line);
	    $sp1=$s[0];
	    $sp2=$s[1];
	}
	else{
	    if($sp1 eq "NA" || $sp2 eq "NA"){
		print "Weird!!! we don't know what the species are: ".$line."\n";
		exit(1);
	    }
	    
	    my @s=split("\t",$line);

	    my $codon1=$s[0];
	    my $codon2=$s[2];
	    my $count=$s[4]+0;
	    
	    if(exists $matrices->{$sp1}){
		if(exists $matrices->{$sp1}{$sp2}){
		    if(exists $matrices->{$sp1}{$sp2}{$codon1}){
			$matrices->{$sp1}{$sp2}{$codon1}{$codon2}=$count;
		    }
		    else{
			$matrices->{$sp1}{$sp2}{$codon1}={$codon2=>$count};
		    }
		}
		else{
		    $matrices->{$sp1}{$sp2}={$codon1=>{$codon2=>$count}};
		}
	    }
	    else{
		$matrices->{$sp1}={$sp2=>{$codon1=>{$codon2=>$count}}};
	    }
	}
	
	$line=<$input>;
    }

    close($input);


    ## check matrices

    my $nbsp=keys %{$matrices};

    foreach my $sp1 (keys %{$matrices}){
	my $nbsp2=keys %{$matrices->{$sp1}};

	print "We have matrices for ".$sp1." and ".$nbsp2." other species.\n";
	
	foreach my $sp2 (keys %{$matrices->{$sp1}}){
	    my $nbcodons=keys %{$matrices->{$sp1}{$sp2}};
	    
	    if($nbcodons<63){
		print "Weird!! we have ".$nbcodons." codons for ".$sp1." and ".$sp2."\n";
		exit(1);
	    }
	    
	    foreach my $c1 (keys %{$matrices->{$sp1}{$sp2}}){
		my $nbcodons2=keys %{$matrices->{$sp1}{$sp2}{$c1}};
		
		if($nbcodons2<63){
		    print "Weird!! we have ".$nbcodons2." codons for ".$sp1." and ".$sp2." and ".$c1."\n";
		    exit(1);
		}
	    }
	}
    }
    
    print "We have matrices for ".$nbsp." species.\n";
}

##############################################################

sub normalizeMatrices{
    my $rawmat=$_[0];
    my $pseudofreq=$_[1];
    my $normmat=$_[2];

    foreach my $sp1 (keys %{$rawmat}){
	$normmat->{$sp1}={};
	
	foreach my $sp2 (keys %{$rawmat->{$sp1}}){
	    $normmat->{$sp1}{$sp2}={};
	    
	    my @codons=keys %{$rawmat->{$sp1}{$sp2}};
	    
	    my $minfreq=1;

	    foreach my $codon1 (@codons){
		$normmat->{$sp1}{$sp2}{$codon1}={};
		
		my $sumfreq=0;
		
		foreach my $codon2 (@codons){
		    $sumfreq+=$rawmat->{$sp1}{$sp2}{$codon1}{$codon2};
		}

		foreach my $codon2 (@codons){
		    my $inifreq=$rawmat->{$sp1}{$sp2}{$codon1}{$codon2};
		    
		    if($sumfreq!=0){
			my $normfreq=($inifreq+0.0)/($sumfreq+0.0);
			$normmat->{$sp1}{$sp2}{$codon1}{$codon2}=$normfreq;
			
			if($normfreq<$minfreq && $normfreq>0){
			    $minfreq=$normfreq;
			}
		    }
		    else{
			$normmat->{$sp1}{$sp2}{$codon1}{$codon2}=$pseudofreq;
		    }
		}
	    }
	}
    }
}

##############################################################

sub removeGapsRef{
    my $aln=$_[0];
    my $refsp=$_[1];
    my $newaln=$_[2];
    
    if(exists $aln->{$refsp}){
	my $size=length $aln->{$refsp}{"sequence"};
       
	
	foreach my $sp (keys %{$aln}){
	    $newaln->{$sp}={"chr"=>$aln->{$sp}{"chr"}, "start"=>$aln->{$sp}{"start"}, "strand"=>$aln->{$sp}{"strand"},"sequence"=>""};
	}
	
	for(my $i=0; $i<$size; $i++){
	    my $refbase=substr $aln->{$refsp}{"sequence"}, $i, 1;
	    
	    if($refbase ne "-"){
		foreach my $sp (keys %{$aln}){
		    my $base=substr $aln->{$sp}{"sequence"}, $i, 1;
		    $newaln->{$sp}{"sequence"}.=$base;
		}
	    }
	}
    }
}

##############################################################

sub computeCSFCodons{
    my $aln=$_[0]; ## there are no gaps in the reference species
    my $refsp=$_[1];
    my $matricescoding=$_[2];
    my $matricesnoncoding=$_[3];
    my $geneticcode=$_[4];
    my $mininformants=$_[5];
    my $penalty=$_[6];
    my $csfscores=$_[7]; ## we keep the score and the number of informative codons for each window
    my $forbidden=$_[8];
    
    my $lengthaln=length $aln->{$refsp}{"sequence"};
    my $chr=$aln->{$refsp}{"chr"};
    my $strand=$aln->{$refsp}{"strand"};
    my $start=$aln->{$refsp}{"start"};
    
    my $increment=1;
    my $pos=$start-1;
    
    if($strand eq "-"){
	$increment=-1;
	$pos=$start+1;
    }

    my @species=keys %{$aln};

    for(my $i=0; $i<$lengthaln; $i++){
	$pos=$pos+$increment;

	my $key=$chr.",".($pos+1); ## middle pos
	
	if(exists $csfscores->{$key}){
	    ## multiple alignments
	    $forbidden->{$key}=1;
	}
	else{
	    my $codonref1=substr $aln->{$refsp}{"sequence"}, $i, 3;
	    $codonref1=uc $codonref1;
	    my $codonref2=reverseComplement($codonref1);
	    
	    my @scores1;
	    my @scores2;
	    
	    if(exists $geneticcode->{$codonref1} &&  exists $geneticcode->{$codonref2}){
		foreach my $sp (@species){
		    if($sp ne $refsp){
			my $codoninf1=substr $aln->{$sp}{"sequence"}, $i, 3;
			$codoninf1=uc $codoninf1;
			my $codoninf2=reverseComplement($codoninf1);
			
			if((exists $geneticcode->{$codoninf1}) && (exists $geneticcode->{$codoninf2})){
			    
			    ## first strand

			    if($codonref1 ne $codoninf1 &&  (exists $matricescoding->{$refsp}{$sp}{$codonref1}{$codoninf1})){
			
				if($matricesnoncoding->{$refsp}{$sp}{$codonref1}{$codoninf1}>0){
				    my $ratio=($matricescoding->{$refsp}{$sp}{$codonref1}{$codoninf1}+0.0)/($matricesnoncoding->{$refsp}{$sp}{$codonref1}{$codoninf1}+0.0);
				    my $thisscore;
				    
				    if($ratio>0){
					$thisscore=log($ratio);
				    }
				    else{
					$thisscore=log($penalty);
				    }
				    
				    push(@scores1,$thisscore);
				}
			    }
			   
			    
			    ## opposite strand
			    
			    if($codonref2 ne $codoninf2  && (exists $matricescoding->{$refsp}{$sp}{$codonref2}{$codoninf2})){
				if($matricesnoncoding->{$refsp}{$sp}{$codonref2}{$codoninf2}>0){
				    my $ratio=($matricescoding->{$refsp}{$sp}{$codonref2}{$codoninf2}+0.0)/($matricesnoncoding->{$refsp}{$sp}{$codonref2}{$codoninf2}+0.0);
				    my $thisscore;
				    
				    if($ratio>0){
					$thisscore=log($ratio);
				    }
				    else{
					$thisscore=log($penalty);
				    }
				
				    push(@scores2,$thisscore);
				}
			    }
			}
		    }
		}
	    }
	    
	    
	    my $nbscores1=@scores1;
	    my $nbscores2=@scores2;

	    if($nbscores1>=$mininformants){
		my $medianscore1=computeMedian(\@scores1);
	    
		if(!(exists $csfscores->{$key})){
		    $csfscores->{$key}={};	
		}
		
		if($strand eq "+"){
		    $csfscores->{$key}{"fwd"}=$medianscore1;
		}
		else{
		    if($strand eq "-"){
			$csfscores->{$key}{"rev"}=$medianscore1;
		    }
		    else{
			print "Weird strand!!! ".$strand."\n";
			exit(1);
		    }
		}
	    }
	    
	    if($nbscores2>=$mininformants){
		my $medianscore2=computeMedian(\@scores2);
		
		if(!(exists $csfscores->{$key})){
		    $csfscores->{$key}={};	
		}
		
		if($strand eq "-"){
		    $csfscores->{$key}{"fwd"}=$medianscore2;
		}
		else{
		    if($strand eq "+"){
			$csfscores->{$key}{"rev"}=$medianscore2;
		    }
		    else{
			print "Weird strand!!! ".$strand."\n";
			exit(1);
		    }
		}
	    }
	}
    }
}

##############################################################

sub computeMedian{
    my $array=$_[0];

    my @sorted=sort {$a<=>$b} @{$array};
    my $nb=@{$array};

    my $median="NA";

    if($nb%2==1){
	$median=$sorted[int($nb/2)];
    }
    else{
	$median=($sorted[(int($nb/2)-1)]+$sorted[int($nb/2)])/2;
    }

    return $median;
}

##############################################################

sub readSelectedInformants{
    my $pathin=$_[0];
    my $ids=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split(" ",$line);
	my $sp=$s[0];

	$ids->{$sp}=1;
	
	$line=<$input>;
    }
    
    close($input);
}

#################################################################

sub reverseComplement{
    my $sequence=$_[0];
    
    my $rev=reverse $sequence;

    $rev=~s/A/X/g;
    $rev=~s/C/Y/g;
    $rev=~s/G/Z/g;
    $rev=~s/T/W/g;

    $rev=~s/X/T/g;
    $rev=~s/Y/G/g;
    $rev=~s/Z/C/g;
    $rev=~s/W/A/g;

    return $rev;
}

##############################################################

sub readSymmetricSubstitutions{
    my $pathin=$_[0];
    my $symm=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);
	my $codon1=$s[0];
	my $codon2=$s[1];
	
	if(exists $symm->{$codon1}){
	    $symm->{$codon1}{$codon2}=1;
	}
	else{
	    $symm->{$codon1}={$codon2=>1};
	}
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub removeSymmetricSubstitutions{
    my $matrices=$_[0];
    my $forbiddenpairs=$_[1];

    my $nbforbidden=keys %{$forbiddenpairs};

    print "There are ".($nbforbidden/2)." forbidden pairs of codons.\n";
     
    foreach my $sp1 (keys %{$matrices}){
	foreach my $sp2 (keys %{$matrices->{$sp1}}){
	    foreach my $codon1 (keys %{$forbiddenpairs}){
		foreach my $codon2  (keys %{$forbiddenpairs->{$codon1}}){
		    delete $matrices->{$sp1}{$sp2}{$codon1}{$codon2};
		}
		
		my $nbleft=keys %{$matrices->{$sp1}{$sp2}{$codon1}};
		
		if($nbleft==0){
		    print "Removing all substitutions for ".$codon1."!!!\n";
		    
		    delete $matrices->{$sp1}{$sp2}{$codon1}; 
		}
	    }
	}
    }
}

##############################################################

sub makeBlocks{
    my $hashcoords=$_[0];
    my $startblocks=$_[1];
    my $endblocks=$_[2];

    my @uniquestart=keys %{$hashcoords};
    my @sortedstart=sort {$a<=>$b} @uniquestart;

    my $nb=@sortedstart;
    
    my $currentstart=$sortedstart[0];
    my $currentend=$hashcoords->{$sortedstart[0]};

    for(my $i=1; $i<$nb; $i++){
	my $thisstart=$sortedstart[$i];
	my $thisend=$hashcoords->{$sortedstart[$i]};

	if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
	    if($thisend>$currentend){
		$currentend=$thisend;
	    }
	}
	else{
	    push(@{$startblocks}, $currentstart);
	    push(@{$endblocks}, $currentend);
	    
	    $currentstart=$thisstart;
	    $currentend=$thisend;
	}
    }
    
    ## don't forget last block
    push(@{$startblocks}, $currentstart);
    push(@{$endblocks}, $currentend);
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes CSF scores on multiple genome alignments.\n";
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
$parameters{"refsp"}="NA";

$parameters{"pathSelectedInformants"}="NA";
$parameters{"pathAlignment"}="NA";

$parameters{"pathCodingMatrices"}="NA";
$parameters{"pathNoncodingMatrices"}="NA";
$parameters{"pathGeneticCode"}="NA";

$parameters{"windowSize"}="NA";
$parameters{"minInformantSpecies"}="NA";
$parameters{"minInformativeCodons"}="NA";
$parameters{"penalty"}=0;
$parameters{"pseudofreq"}=0;

$parameters{"pathSymmetricSubstitutions"}="NA";
$parameters{"includeSymmetric"}="NA";


$parameters{"start"}=0;
$parameters{"end"}=0;

$parameters{"pathOutputPositiveScores"}="NA";
$parameters{"pathOutputCoveredRegions"}="NA";

my @defaultpars=("refsp", "pathSelectedInformants", "pathAlignment","pathCodingMatrices","pathNoncodingMatrices","pathGeneticCode","windowSize","minInformantSpecies","minInformativeCodons","penalty","pseudofreq", "pathSymmetricSubstitutions", "includeSymmetric", "start", "end", "pathOutputPositiveScores", "pathOutputCoveredRegions");

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

print "Reading genetic code...\n";

my %geneticcode;
readGeneticCode($parameters{"pathGeneticCode"},\%geneticcode);
my $nbcodons=keys %geneticcode;
print "There are ".$nbcodons." codons in the genetic code.\n";

print "Done.\n\n";

##############################################################

print "Reading matrices...\n";

my %codingraw;
readMatrices($parameters{"pathCodingMatrices"},\%codingraw);

my %noncodingraw;
readMatrices($parameters{"pathNoncodingMatrices"},\%noncodingraw);

 print "Done.\n\n";


my $includeSymmetric=$parameters{"includeSymmetric"};

if($includeSymmetric eq "no"){

    print "Reading symmetric substitutions...\n";

    my %symm;
    readSymmetricSubstitutions($parameters{"pathSymmetricSubstitutions"},\%symm);

    print "Done.\n\n";
    
    print "Removing symmetric synonymous substitutions from the matrices...\n";
    
    removeSymmetricSubstitutions(\%codingraw,\%symm);
    
    removeSymmetricSubstitutions(\%noncodingraw,\%symm);
  
    print "Done.\n\n";
}
else{
    print "We keep symmetric substitutions.\n";
}


print "Normalizing matrices...\n";
my $pseudofreq=$parameters{"pseudofreq"};

print "Pseudo frequency: ".$pseudofreq."\n";
my %codingnorm;
normalizeMatrices(\%codingraw, $pseudofreq, \%codingnorm);

my %noncodingnorm;
normalizeMatrices(\%noncodingraw,$pseudofreq,\%noncodingnorm);

print "Done.\n\n";

##############################################################

print "Reading selected informants...\n";

my %informants;

readSelectedInformants($parameters{"pathSelectedInformants"},\%informants);

my $nbinf=keys %informants;

print "Kept ".$nbinf." informant species: ".join(", ", keys %informants)."\n";

print "Done.\n";

##############################################################

print "Reading alignments and computing CSF scores for each codon...\n";

my $start=$parameters{"start"}+0;
my $end=$parameters{"end"}+0;

print "start ".$start." end ".$end."\n";

my $refsp=$parameters{"refsp"};

print "Reference species: ".$refsp."\n";

my $mininformants=$parameters{"minInformantSpecies"}+0;
my %codoncsf;
my %forbidden;
my %chrsizesref;

my $penalty=$parameters{"penalty"}+0.0;

readAlignmentsCSF($parameters{"pathAlignment"}, $refsp, \%informants, \%codingnorm, \%noncodingnorm, \%geneticcode, $penalty, $mininformants, $start, $end, \%codoncsf, \%forbidden, \%chrsizesref);
 
my $nbcodons=keys %codoncsf;
my $nbforbidden=keys %forbidden;

foreach my $key (keys %forbidden){
    delete $codoncsf{$key};
}

print "We computed the score for ".$nbcodons." codons, including ".$nbforbidden." forbidden positions (multiple alignments).\n";

print "Done.\n";

##############################################################

print "Computing score for windows and writing output...\n";
my $windowSize=$parameters{"windowSize"}+0;
my $mininfocodons=$parameters{"minInformativeCodons"}+0;

print "window size ".$windowSize." min info codons ".$mininfocodons."\n";

open(my $outputcoverage, ">".$parameters{"pathOutputCoveredRegions"});
print $outputcoverage "Chr\tStart\tEnd\tStrand\n";

open(my $output, ">".$parameters{"pathOutputPositiveScores"});

print $output "Chr\tStart\tEnd\tStrand\tSumScore\tNbInfoCodons\n";

foreach my $chr (keys %chrsizesref){
    my $size=$chrsizesref{$chr};
    my $limit=min(($end+$windowSize),$size);

    my %hashcoveredplus;
    my %hashcoveredminus;
    
    for(my $i=$start; $i<$limit; $i++){
	
	my $startwin=$i+1;
	my $endwin=$i+$windowSize;

	my @scoresplus;
	my @scoresminus;

	for(my $pos=$startwin+1; $pos<$endwin; $pos+=3){
	    my $key=$chr.",".$pos;

	    if(exists $codoncsf{$key}){
		if(exists $codoncsf{$key}{"fwd"}){
		    push(@scoresplus, $codoncsf{$key}{"fwd"});
		}

		if(exists $codoncsf{$key}{"rev"}){
		    push(@scoresminus, $codoncsf{$key}{"rev"});
		}
	    }
	}
	
	my $nbinfoplus=@scoresplus;
	my $nbinfominus=@scoresminus;
	
	if($nbinfoplus>=$mininfocodons){
	    my $sumplus=sum @scoresplus;
	    my $sp=sprintf("%.2f",$sumplus);

	    if($sumplus>0){
		print $output $chr."\t".$startwin."\t".$endwin."\t+\t".$sp."\t".$nbinfoplus."\n";
	    }

	    if(exists $hashcoveredplus{$startwin}){
		print "Weird! already saw start window ".$startwin." for fwd strand.\n";
		exit(1);
	    }

	    $hashcoveredplus{$startwin}=$endwin;
	}

	if($nbinfominus>=$mininfocodons){
	    my $summinus=sum @scoresminus;
	  
	    my $sm=sprintf("%.2f",$summinus);
	    
	    if($summinus>0){
		print $output $chr."\t".$startwin."\t".$endwin."\t-\t".$sm."\t".$nbinfominus."\n";
	    }

	    if(exists $hashcoveredminus{$startwin}){
		print "Weird! already saw start window ".$startwin." for reverse strand.\n";
		exit(1);
	    }

	    $hashcoveredminus{$startwin}=$endwin;
	}
    }

    print "Making output covered regions for ".$chr."...\n";

    ## make blocks 
    my @startplus;
    my @endplus;

    print "Forward strand\n";

    makeBlocks(\%hashcoveredplus, \@startplus, \@endplus);

    for(my $i=0; $i<@startplus; $i++){
	print $outputcoverage $chr."\t".$startplus[$i]."\t".$endplus[$i]."\t+\n";
    }

    print "Reverse strand\n";

    my @startminus;
    my @endminus;

    makeBlocks(\%hashcoveredminus, \@startminus, \@endminus);
    
    for(my $i=0; $i<@startminus; $i++){
	print $outputcoverage $chr."\t".$startminus[$i]."\t".$endminus[$i]."\t-\n";
    }
    
    print "Done.\n";
    
}

close($output);
close($outputcoverage);

print "Done.\n";

##############################################################
