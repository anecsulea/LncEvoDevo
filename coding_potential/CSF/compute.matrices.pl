use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX;
use strict;

##############################################################

sub readRegionCoordinates{
    my $pathin=$_[0];
    my $coords=$_[1];

    open(my $input, $pathin);
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    my %header;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    my %unordered;

    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $gene=$s[$header{"Gene"}];
	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	my $frame=$s[$header{"Frame"}]+0;

	if($frame!=0){
	    if($strand eq "1"){
		$start=$start+(3-$frame);
	    }
	    else{
		if($strand eq "-1"){
		    $start=$start+$frame;
		}
		else{
		    print "Weird strand for ".$gene."!!!\n";
		    exit(1);
		}
	    }
	}
	
	my $length=$end-$start+1;
	my $reste=$length%3;

	$end=$end-$reste;

	my $efflen=$end-$start+1;

	if($efflen>=30){ 
	    ## at least 10 codons
	    
	    if(exists $unordered{$chr}){
		if(exists $unordered{$chr}{$start}){
		    push(@{$unordered{$chr}{$start}{"end"}}, $end);
		    push(@{$unordered{$chr}{$start}{"strand"}}, $strand);
		    push(@{$unordered{$chr}{$start}{"gene"}}, $gene);
		}
		else{
		    $unordered{$chr}{$start}={"gene"=>[$gene], "end"=>[$end], "strand"=>[$strand]};
		}
	    }
	    else{
		$unordered{$chr}={$start=>{"gene"=>[$gene], "end"=>[$end], "strand"=>[$strand]}};
	    }
	}
	    
	$line=<$input>;
    }
    
    close($input);

    my $nbreg=0;

    foreach my $chr (keys %unordered){
	$coords->{$chr}={"start"=>[], "end"=>[], "strand"=>[], "gene"=>[]};

	my @uniquestart=keys %{$unordered{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nbe=@{$unordered{$chr}{$start}{"end"}};

	    $nbreg+=$nbe;

	    for(my $i=0; $i<$nbe; $i++){
		push(@{$coords->{$chr}{"start"}}, $start);
		push(@{$coords->{$chr}{"end"}}, ${$unordered{$chr}{$start}{"end"}}[$i]);
		push(@{$coords->{$chr}{"strand"}}, ${$unordered{$chr}{$start}{"strand"}}[$i]);
		push(@{$coords->{$chr}{"gene"}}, ${$unordered{$chr}{$start}{"gene"}}[$i]);
	    }
	}
    }

    print "Kept ".$nbreg." regions.\n";
}

##############################################################

sub constructRegionIndex{
    my $orderedregions=$_[0];
    my $windowsize=$_[1];
    my $regionindex=$_[2];

    $regionindex->{"window"}=$windowsize+0;

    foreach my $chr (keys %{$orderedregions}){
	$regionindex->{$chr}={};

	my $nbregions=@{$orderedregions->{$chr}{"start"}};

	$regionindex->{$chr}{0}=0;
	
	for(my $i=0; $i<$nbregions; $i++){
	    my $end=${$orderedregions->{$chr}{"end"}}[$i];

	    my $pos=$windowsize*(ceil($end/$windowsize));
	    
	    $regionindex->{$chr}{$pos}=$i; ## last index i for which end <= $pos
	}
    }
}

##############################################################

sub readAlignments{
    my $pathin=$_[0];
    my $refsp=$_[1];
    my $regioncoords=$_[2];
    my $regionindex=$_[3];
    my $minlen=$_[4];
    my $minspecies=$_[5];
    my $matrices=$_[6];

    ## start reading the alignments
    ## the alignments are not necessarily ordered ! but we have an index in the regions to know where to look, efficiently
    
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
 
    my $inaln=0;
    my %aln;
    my $nbdone=0;

    my $currentstart="NA";
   
    while($line){

	chomp $line;
	$firstchar=substr $line,0,1;
	
	if($firstchar eq "a"){
	    if($inaln==1){
		## we've already seen an alignment, we check if it overlaps with a region

		if(exists $aln{$refsp}){
		    
		    my $chr=$aln{$refsp}{"chr"};
		    
		    if(exists $regioncoords->{$chr}){
			my $nbregions=@{$regioncoords->{$chr}{"start"}};
			
			my $refstart=$aln{$refsp}{"start"};
			my $reflen=$aln{$refsp}{"length"};
			my $refstrand=$aln{$refsp}{"strand"};
			my $refend=$refstart+$reflen-1;
			my $nbsp=keys %aln;
			
			if($refstrand eq "-"){
			    $refend=$aln{$refsp}{"start"};
			    $refstart=$refend-$reflen+1;
			    ## this way refstart < refend
			}
			
			## check if this alignment overlaps with any region
			
			my $floorindex=floor($refstart/$regionindex->{"window"})*$regionindex->{"window"};
			
			if(exists $regionindex->{$chr}{$floorindex}){
			    
			    my $i=$regionindex->{$chr}{$floorindex}; ## this position should be before refstart, we don't miss anything by not looking before
			    
			    while($i<$nbregions && ${$regioncoords->{$chr}{"end"}}[$i]<$refstart){
				$i++;
			    }
			    
			    ## if there are enough species and we have minimum length, we check for overlapping regions
			    
			    if($nbsp>=$minspecies && $reflen>=$minlen){
				my %ovregions;
				$ovregions{"start"}=[];
				$ovregions{"end"}=[];
				$ovregions{"strand"}=[];
				
				while($i<$nbregions && ${$regioncoords->{$chr}{"start"}}[$i]<=$refend){
				    my $M=max($refstart, ${$regioncoords->{$chr}{"start"}}[$i]);
				    my $m=min($refend, ${$regioncoords->{$chr}{"end"}}[$i]);
				    
				    if($M<=$m){
					push(@{$ovregions{"start"}}, ${$regioncoords->{$chr}{"start"}}[$i]);
					push(@{$ovregions{"end"}}, ${$regioncoords->{$chr}{"end"}}[$i]);
					push(@{$ovregions{"strand"}}, ${$regioncoords->{$chr}{"strand"}}[$i]);
				    }
				    
				    $i++;
				}
				
				my $nbov=@{$ovregions{"start"}};
				
				if($nbov>0){
				    ## print "Analyzing aln ".$refsp." ".$chr." ".$refstart." ".$refstrand." length ".$reflen." ".$nbov." overlapping regions \n";
				    
				    computeMatrices(\%aln, $refsp, \%ovregions, $matrices);
				}
			    }
			}
		    }
		}
		
	    }

	    ## reinitialize alignment
	    
	    %aln=();
	    $inaln=1;
	    $nbdone++;

	    if($nbdone%10000==0){
		print $nbdone." alignments done.\n";
	    }
	}
	
	if($firstchar eq "s"){
	    my @s=split(" ",$line);
	    my $spinfo=$s[1];
	    my @t=split("\\.", $spinfo);
	    my $sp=$t[0];
	    my $chr=$t[1];
	    
	    my $start=$s[2]+1; ## is now 1-based 
	    my $ungappedlen=$s[3];
	    my $strand=$s[4];
	    my $sequence=$s[6];
	    my $lenchr=$s[5]+0;

	    if($strand eq "-"){
		$start=$lenchr-$s[2];
	    }

	    $aln{$sp}={"chr"=>$chr, "start"=>$start, "length"=>$ungappedlen, "strand"=>$strand,"sequence"=>$sequence};
	}
		
	$line=<$input>;
    }

    close($input);
 }

##############################################################

sub computeMatrices{
    my $aln=$_[0];
    my $refsp=$_[1];
    
    my $regioncoords=$_[2];
    my $matrices=$_[3];    
       
    ## first, read the regions and make a hash of codon positions
    
    my %hashcodons;
    my %dualstrand;
    
    my $nbreg=@{$regioncoords->{"start"}};
    
    for(my $i=0;  $i<$nbreg; $i++){
	my $regstart=${$regioncoords->{"start"}}[$i];
	my $regend=${$regioncoords->{"end"}}[$i];
	my $regstrand=${$regioncoords->{"strand"}}[$i];
	
	my $len=$regend-$regstart+1;
	my $mod3=$len%3;

	if($mod3!=0){
	    print "Region length is not multiple of 3! ".$regstart." ".$regend."\n";
	    exit(1);
	}	
	
	for(my $j=$regstart; $j<=($regend-2); $j+=3){
	    my $key=$j."-".($j+2);
	    if(exists $hashcodons{$key}){
		if($hashcodons{$key} ne $regstrand){
		    $dualstrand{$key}=1;
		}
	    }
	    else{
		$hashcodons{$key}=$regstrand;
	    }
	}
    }

    my $nbdual=keys %dualstrand;
    if($nbdual>0){
	print "Found coding regions on both strands ".join(";", @{$regioncoords->{"start"}})." ".join(";", @{$regioncoords->{"end"}})."\n";

	foreach my $j (keys %dualstrand){
	    delete $hashcodons{$j};
	}
    }

    ### now go through the alignment and add to matrices

    my $refpos=$aln->{$refsp}{"start"}-1;
    my $increment=1;

    my $refstrand=$aln->{$refsp}{"strand"};

    if($refstrand eq "-"){
	$refpos=$aln->{$refsp}{"start"}+1;
	$increment=-1;
    } else{
	if($refstrand ne "+"){
	    print "Weird strand for reference species: ".$refstrand.".\n";
	    exit(1);
	}
    }

    my $lengthseq=length $aln->{$refsp}{"sequence"};

    for(my $i=0; $i<($lengthseq-2); $i++){
	my $base1ref=uc (substr $aln->{$refsp}{"sequence"}, $i, 1);
	
	if($base1ref ne "-"){
	    $refpos+=$increment;

	    my $keypos=$refpos."-".($refpos+2);

	    if($refstrand eq "-"){ ## coordinates go in the opposite direction
		$keypos=($refpos-2)."-".$refpos;
	    }

	    if(exists $hashcodons{$keypos}){
		my $strandcodon=$hashcodons{$keypos};

		if($strandcodon ne "1" && $strandcodon ne "-1"){
		    print "Weird codon strand for ".$keypos."!!\n";
		    exit(1);
		}

		## check if there are gaps afterwards
		
		my $base2ref=uc (substr $aln->{$refsp}{"sequence"}, ($i+1), 1);
		my $base3ref=uc (substr $aln->{$refsp}{"sequence"}, ($i+2), 1);
		
		if($base2ref ne "-" && $base3ref ne "-" && $base1ref ne "N" && $base2ref ne "N"  && $base3ref ne "N"){
		    my $codonref=$base1ref.$base2ref.$base3ref;

		    if(($strandcodon eq "-1" && $refstrand eq "+") || ($strandcodon eq "1" && $refstrand eq "-")){
			$codonref=reverseComplement($codonref); ## we do this just once for the reference species
		    }
		    
		    foreach my $othersp (keys %{$aln}){
			if($othersp ne $refsp){
			    my $codonother=uc (substr $aln->{$othersp}{"sequence"}, $i, 3);
			    my $countA = ($codonother =~ tr/A//);
			    my $countC = ($codonother =~ tr/C//);
			    my $countG = ($codonother =~ tr/G//);
			    my $countT = ($codonother =~ tr/T//);
			    
			    my $tot=$countA+$countC+$countG+$countT;
			    
			    if($tot==3){
				## there are no gaps in either of the species

				if(($strandcodon eq "-1" && $refstrand eq "+") || ($strandcodon eq "1" && $refstrand eq "-")){
				    $codonother=reverseComplement($codonother);
				}
				
				if(exists $matrices->{$othersp}){
				    if(exists $matrices->{$othersp}{$codonref}){
					if(exists $matrices->{$othersp}{$codonref}{$codonother}){
					    $matrices->{$othersp}{$codonref}{$codonother}++;
					}
					else{
					    $matrices->{$othersp}{$codonref}{$codonother}=1;
					}
				    }
				    else{
					$matrices->{$othersp}{$codonref}={$codonother=>1};
				    }
				}
				else{
				    $matrices->{$othersp}={$codonref=>{$codonother=>1}};
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

##############################################################

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

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes CSF matrices.\n";
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

$parameters{"pathAlignment"}="NA";
$parameters{"refSpecies"}="NA";
$parameters{"minAlignedLength"}=0;
$parameters{"minAlignedSpecies"}=0;
$parameters{"pathRegionCoordinates"}="NA";
$parameters{"pathGeneticCode"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathAlignment", "refSpecies", "minAlignedLength", "minAlignedSpecies", "pathRegionCoordinates", "pathGeneticCode", "pathOutput");


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

print "Reading region coordinates...\n";

my %regioncoords;
readRegionCoordinates($parameters{"pathRegionCoordinates"},\%regioncoords);

print "Done.\n";

print "Constructing region index, 100kb window...\n";
my %regionindex;

constructRegionIndex(\%regioncoords, 100000, \%regionindex);

print "Done.\n";

##############################################################

print "Reading genetic code...\n";

my %geneticcode;
readGeneticCode($parameters{"pathGeneticCode"},\%geneticcode);
my $nbcodons=keys %geneticcode;
print "There are ".$nbcodons." codons in the genetic code.\n";

print "Done.\n\n";

##############################################################

print "Reading alignments and computing matrices...\n";

my $refsp=$parameters{"refSpecies"};
my $minalnlen=$parameters{"minAlignedLength"}+0;
my $minalnsp=$parameters{"minAlignedSpecies"}+0;

print "Reference species: ".$refsp."\n";
print "Minimum aligned length: ".$minalnlen."bp.\n";
print "Minimum aligned species: ".$minalnsp.".\n";

my %matrices;

readAlignments($parameters{"pathAlignment"},  $refsp, \%regioncoords, \%regionindex, $minalnlen, $minalnsp,  \%matrices);

print "Done.\n\n";

################################################################


print "Writing output...\n";

my @codons=keys %geneticcode;
    
open(my $output,">".$parameters{"pathOutput"});

foreach my $sp2 (keys %matrices){
    print $output "# ".$refsp." ".$sp2."\n";
    
    foreach my $codon1 (@codons){
	foreach my $codon2 (@codons){
	    if(exists $matrices{$sp2}{$codon1}{$codon2}){
		print $output $codon1."\t".$geneticcode{$codon1}."\t".$codon2."\t".$geneticcode{$codon2}."\t".$matrices{$sp2}{$codon1}{$codon2}."\n";
	    }
	    else{
		print $output $codon1."\t".$geneticcode{$codon1}."\t".$codon2."\t".$geneticcode{$codon2}."\t"."\t0\n";
	    }
	}
    }
}

close($output);

print "Done.\n\n";

###############################################################
