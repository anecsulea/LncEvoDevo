use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $exons=$_[1];
    my $genestx=$_[2];
    my $txexons=$_[3];

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
		} else{
		    print "Weird strand!\n";
		    print $line."\n";
		    exit(1);
		} 
	    }

	    my $info=$s[8];
	    my @infoarray=split(";", $info);

	    my $geneid=findInfo("gene_id", \@infoarray);
	    my $txid=findInfo("transcript_id", \@infoarray);
	    my $exonid=findInfo("exon_id", \@infoarray);
	    
	    if($exonid eq "NA"){
		$exonid=$chr.",".$start.",".$end.",".$strand;
	    }

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

############################################################################

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

############################################################################

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
	$refordered->{$chr}={"start"=>[], "end"=>[], "strand"=>[], "id"=>[]};
	
	my @uniquestart=keys %{$hashstart{$chr}};
	
	my @sortedstart = sort {$a <=> $b} @uniquestart;
	
	foreach my $b (@sortedstart){
	    
	    my $nbblocks=@{$hashstart{$chr}{$b}{"end"}};
	    
	    for(my $i=0;$i<$nbblocks;$i++){
		push(@{$refordered->{$chr}{"start"}},$b);
		push(@{$refordered->{$chr}{"strand"}},${$hashstart{$chr}{$b}{"strand"}}[$i]);
		push(@{$refordered->{$chr}{"end"}},${$hashstart{$chr}{$b}{"end"}}[$i]);
		push(@{$refordered->{$chr}{"id"}},${$hashstart{$chr}{$b}{"id"}}[$i]);
	    }
	}
    }
}

######################################################################

sub readRepeatMasker{
    my $paths=$_[0];
    my $repeatmasker=$_[1];
    
    my %hashrepeats;

    foreach my $pathin (@{$paths}){
	print "Reading from ".$pathin."\n";

	my @s=split("\\.", $pathin);
	my $ext=$s[-1];
	my $input;
	
	if($ext eq "gz"){
	    open($input,"zcat $pathin |");
	} else{
	    open($input,$pathin);
	}
	
	
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
	    my @s=split(" ",$line);
	    
	    my $chr=$s[5];
	    
	    my $prefix=substr $chr,0,3;
	    
	    if($prefix eq "chr"){
		$chr=substr $chr,3;
	    }
	    
	    
	    my $start=$s[6]+1; ## coordinates start at 1 now
	    my $end=$s[7]+0;
	    my $strand=$s[9];

	    if($strand eq "+"){
		$strand="1";
	    }
	    else{
		$strand="-1";
	    }
	    
	    if(exists $hashrepeats{$chr}){
		if(exists $hashrepeats{$chr}{$start}){
		    push(@{$hashrepeats{$chr}{$start}{"end"}},$end);
		    push(@{$hashrepeats{$chr}{$start}{"strand"}},$strand);
		}
		else{
		    $hashrepeats{$chr}{$start}={"end"=>[$end], "strand"=>[$strand]};
		}
	    }
	    else{
		$hashrepeats{$chr}={$start=>{"end"=>[$end],"strand"=>[$strand]}};
	    }
	    
	    
	    $line=<$input>;
	}
	
	close($input);
    }
    
    ### now order repeats
    
    foreach my $chr (keys %hashrepeats){
	$repeatmasker->{$chr}={"start"=>[],"end"=>[],"strand"=>[]};
	
	my @uniquestart=keys %{$hashrepeats{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;
	
	foreach my $start (@sortedstart){
	    my $nb=@{$hashrepeats{$chr}{$start}{"end"}};
	    for(my $i=0; $i<$nb; $i++){
		
		push(@{$repeatmasker->{$chr}{"start"}},$start);
		push(@{$repeatmasker->{$chr}{"end"}},${$hashrepeats{$chr}{$start}{"end"}}[$i]);
		push(@{$repeatmasker->{$chr}{"strand"}},${$hashrepeats{$chr}{$start}{"strand"}}[$i]);
	    }
	}
    }
}

##############################################################

sub readRetrogenes{
    my $paths=$_[0];
    my $retrogenes=$_[1];

    my %hashretro;
    
    foreach my $pathin (@{$paths}){
	print "Reading from ".$pathin."\n";

	my @s=split("\\.", $pathin);
	my $ext=$s[-1];
	
	my $input;
	
	if($ext eq "gz"){
	    open($input, "zcat $pathin |");
	} else{
	    open($input, $pathin);
	}
	
	my $line=<$input>;

	while($line){
	    chomp $line;
	    my @s=split("\t", $line);
	    my $chr=$s[0];
	    my $prefix=substr $chr, 0,3;
	    
	    if($prefix eq "chr"){
		$chr=substr $chr,3;
	    }
	    
	    my $start=$s[1]+0;
	    my $end=$s[2]+0;
	    
	    my $ss=$s[5];
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
	    
	    if(exists $hashretro{$chr}){
		if(exists $hashretro{$chr}{$start}){
		    push(@{$hashretro{$chr}{$start}{"end"}},$end);
		    push(@{$hashretro{$chr}{$start}{"strand"}},$strand);
		}
		else{
		    $hashretro{$chr}{$start}={"end"=>[$end], "strand"=>[$strand]};
		}
	    }
	    else{
		$hashretro{$chr}={$start=>{"end"=>[$end],"strand"=>[$strand]}};
	    }
	    

	    $line=<$input>;
	}
	
	close($input);
    }

    ### now order retrogenes
    
    foreach my $chr (keys %hashretro){
	$retrogenes->{$chr}={"start"=>[],"end"=>[],"strand"=>[]};
	
	my @uniquestart=keys %{$hashretro{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashretro{$chr}{$start}{"end"}};
	    for(my $i=0; $i<$nb; $i++){

		push(@{$retrogenes->{$chr}{"start"}},$start);
		push(@{$retrogenes->{$chr}{"end"}},${$hashretro{$chr}{$start}{"end"}}[$i]);
		push(@{$retrogenes->{$chr}{"strand"}},${$hashretro{$chr}{$start}{"strand"}}[$i]);
	    }
	}
    }
}

##########################################################################

sub extractOverlapExons{
    ## both strands, but not the same gene
    my $exoncoords=$_[0];
    my $exoninfo=$_[1];
    my $overlapcoords=$_[2];

    foreach my $chr (keys %{$exoncoords}){
	my $nbexons=@{$exoncoords->{$chr}{"start"}};
	
	my $firstj=0;

	for(my $i=0; $i<$nbexons; $i++){
	    my $start1=${$exoncoords->{$chr}{"start"}}[$i];
	    my $end1=${$exoncoords->{$chr}{"end"}}[$i];
	    my $id1=${$exoncoords->{$chr}{"id"}}[$i];

	    my $j=$firstj;

	    while($j<$nbexons && ${$exoncoords->{$chr}{"end"}}[$j]<$start1){
		$j++;
	    }
	    
	    $firstj=$j;

	    while($j<$nbexons && ${$exoncoords->{$chr}{"start"}}[$j]<=$end1){
		my $start2=${$exoncoords->{$chr}{"start"}}[$j];
		my $end2=${$exoncoords->{$chr}{"end"}}[$j];
		my $id2=${$exoncoords->{$chr}{"id"}}[$j];

		my $M=max($start1, $start2);
		my $m=min($end1, $end2);
		
		if($M<=$m){
		    my $diffgene=0;
		    
		    foreach my $gene1 (keys %{$exoninfo->{$id1}{"genes"}}){
			foreach my $gene2 (keys %{$exoninfo->{$id2}{"genes"}}){
			    if($gene1 ne $gene2){
				$diffgene=1;
				last;
			    }
			}
		    }
		    
		    if($diffgene==1){
			if(exists $overlapcoords->{$id1}){
			    push(@{$overlapcoords->{$id1}{"start"}}, $M);
			    push(@{$overlapcoords->{$id1}{"end"}}, $m);
			}
			else{
			    $overlapcoords->{$id1}={"start"=>[$M], "end"=>[$m]};
			}
		    }
		}  
		
		$j++;
		
	    }
	}
    }
}

##########################################################################

sub extractOverlapRepeats{
    ## both strands also
    my $exoncoords=$_[0];
    my $repeatcoords=$_[1];
    my $overlapcoords=$_[2];

    foreach my $chr (keys %{$exoncoords}){
	if(exists $repeatcoords->{$chr}){
	    my $nbexons=@{$exoncoords->{$chr}{"start"}};
	    my $nbrepeats=@{$exoncoords->{$chr}{"end"}};
	    
	    my $firstj=0;
	    
	    for(my $i=0; $i<$nbexons; $i++){
		my $start1=${$exoncoords->{$chr}{"start"}}[$i];
		my $end1=${$exoncoords->{$chr}{"end"}}[$i];
		my $id1=${$exoncoords->{$chr}{"id"}}[$i];
		
		my $j=$firstj;
		
		while($j<$nbrepeats && ${$repeatcoords->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}
		
		$firstj=$j;
		
		while($j<$nbrepeats && ${$repeatcoords->{$chr}{"start"}}[$j]<=$end1){
		    my $start2=${$repeatcoords->{$chr}{"start"}}[$j];
		    my $end2=${$repeatcoords->{$chr}{"end"}}[$j];
		    my $id2=${$repeatcoords->{$chr}{"id"}}[$j];
		    
		    my $M=max($start1, $start2);
		    my $m=min($end1, $end2);
		    
		    if($M<=$m){
			if(exists $overlapcoords->{$id1}){
			    push(@{$overlapcoords->{$id1}{"start"}}, $M);
			    push(@{$overlapcoords->{$id1}{"end"}}, $m);
			}
			else{
			    $overlapcoords->{$id1}={"start"=>[$M], "end"=>[$m]};
			}
		    }
		    
		    $j++;
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
		
		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}
	
	push(@{$blocks->{"start"}}, $currentstart);
	push(@{$blocks->{"end"}}, $currentend);
    }
}

##############################################################

sub extractNonOverlappingRegions{
    my $exoninfo=$_[0];
    my $exonoverlap=$_[1];
    my $nonoverlap=$_[2];

    foreach my $id (keys %{$exoninfo}){
	my $start=$exoninfo->{$id}{"start"};
	my $end=$exoninfo->{$id}{"end"};
	
	if(exists $exonoverlap->{$id}){
	    $nonoverlap->{$id}={"start"=>[], "end"=>[]};
	    
	    my %blocksoverlap;
	    makeBlocks($exonoverlap->{$id}, \%blocksoverlap);
	    
	    my $nbov=@{$blocksoverlap{"start"}};

	    my $firststart=${$blocksoverlap{"start"}}[0];
	    if($firststart>$start){
		push(@{$nonoverlap->{$id}{"start"}}, $start);
		push(@{$nonoverlap->{$id}{"end"}}, $firststart-1);
	    }
	    
	    if($nbov>=2){
		for(my $i=0; $i<($nbov-1); $i++){
		    my $thisend=${$blocksoverlap{"end"}}[$i];
		    my $nextstart=${$blocksoverlap{"start"}}[$i+1];

		    if($thisend<($nextstart-1)){
			push(@{$nonoverlap->{$id}{"start"}}, $thisend+1);
			push(@{$nonoverlap->{$id}{"end"}}, $nextstart-1);
		    }
		}
	    }

	    my $lastend=${$blocksoverlap{"end"}}[-1];
	    if($lastend<$end){
		push(@{$nonoverlap->{$id}{"start"}}, $lastend+1);
		push(@{$nonoverlap->{$id}{"end"}}, $end);
	    }
	    
	} else{
	    $nonoverlap->{$id}={"start"=>[$start], "end"=>[$end]};
	}
    }
}
##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts non-overlapping regions for exons, repeats and retrogenes.\n";
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
$parameters{"pathRepeats"}="NA";
$parameters{"pathRetrogenes"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathAnnotGTF", "pathRepeats", "pathRetrogenes", "pathOutput");

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

readGTF($parameters{"pathAnnotGTF"}, \%exons, \%genestx, \%txexons);

my $nbex=keys %exons;
my $nbg=keys %genestx;
my $nbtx=keys %txexons;

print "Found ".$nbex." exons, ".$nbtx." transcripts and ".$nbg." genes.\n";

print "Done.\n";

##############################################################

print "Reading repeats...\n";

my @pathsrep=split(",",$parameters{"pathRepeats"});
my %repeats;

readRepeatMasker(\@pathsrep, \%repeats);

print "Done.\n";


print "Reading repeats...\n";

my @pathsretro=split(",",$parameters{"pathRetrogenes"});
my %retrogenes;

readRetrogenes(\@pathsretro, \%retrogenes);

print "Done.\n";

##############################################################

print "Ordering exons...\n";

my %orderedexons;

orderExons(\%exons, \%orderedexons);

print "Done.\n";

##############################################################

print "Computing overlap...\n";
my %exonoverlap;

print "with other genes: \n";
extractOverlapExons(\%orderedexons, \%exons, \%exonoverlap);
my $nbov=keys %exonoverlap;
print $nbov." exons with overlap.\n";


print "with repeats: \n";
extractOverlapRepeats(\%orderedexons, \%repeats, \%exonoverlap);
my $nbov=keys %exonoverlap;
print $nbov." exons with overlap.\n";

print "with retrogenes: \n";
extractOverlapRepeats(\%orderedexons, \%retrogenes, \%exonoverlap);
my $nbov=keys %exonoverlap;
print $nbov." exons with overlap.\n";

print "Done.\n";

##############################################################

print "Computing non-overlapping regions...\n";
my %nonoverlap;
extractNonOverlappingRegions(\%exons, \%exonoverlap, \%nonoverlap);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tTranscriptID\tChr\tStrand\tExonStart\tExonEnd\tNonoverlapStart\tNonoverlapEnd\n";

foreach my $geneid (keys %genestx){
    foreach my $txid (keys %{$genestx{$geneid}}){
	foreach my $exonid (keys %{$txexons{$txid}}){
	    my $chr=$exons{$exonid}{"chr"};
	    my $start=$exons{$exonid}{"start"};
	    my $end=$exons{$exonid}{"end"};
	    my $strand=$exons{$exonid}{"strand"};
	    
	    if(!exists $nonoverlap{$exonid}){
		print "Weird! cannot find ".$exonid." in non-overlapping regions hash table.\n";
		exit(1);
	    }
	    
	    my $nbnonov=@{$nonoverlap{$exonid}{"start"}};

	    if($nbnonov>0){
		for(my $i=0; $i<$nbnonov; $i++){
		    print $output $geneid."\t".$txid."\t".$chr."\t".$strand."\t".$start."\t".$end."\t".${$nonoverlap{$exonid}{"start"}}[$i]."\t".${$nonoverlap{$exonid}{"end"}}[$i]."\n";
		}
	    } else{
		print $output $geneid."\t".$txid."\t".$chr."\t".$strand."\t".$start."\t".$end."\tNA\tNA\n";
	    }
	}
    }
}

close($output);

print "Done.\n";

##############################################################
