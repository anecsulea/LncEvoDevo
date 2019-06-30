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

sub readRepeatMasker{
    my $pathin=$_[0];
    my $repeats=$_[1];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    my $line=<$input>;
    my %header;
    chomp $line;
    my @s=split("\t", $line);
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $chr=$s[$header{"genoName"}];
	my $prefix=substr $chr, 0,3;

	if($prefix eq "chr"){
	    $chr=substr $chr,3;
	}

	my $start=$s[$header{"genoStart"}]+0;
	my $end=$s[$header{"genoEnd"}]+0;

	my $ss=$s[$header{"strand"}];
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

	$repeats->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	
	$line=<$input>;
    }

    close($input);

}

##############################################################

sub computeOverlap{
    my $orderedexons=$_[0];
    my $repeats=$_[1];
    my $type=$_[2];
    my $overlap=$_[3];

    if($type ne "sense" && $type ne "antisense" && $type ne "any"){
	print "Weird type: ".$type."\n";
	exit(1);
    }
    
    foreach my $chr (keys %{$orderedexons}){
	if(exists $repeats->{$chr}){
	    my $nbex=@{$orderedexons->{$chr}{"start"}};
	    my $nbreg=@{$repeats->{$chr}{"start"}};
	    
	    my $firstreg=0;
	    
	    for(my $i=0; $i<$nbex; $i++){
		my $startex=${$orderedexons->{$chr}{"start"}}[$i];
		my $endex=${$orderedexons->{$chr}{"end"}}[$i];
		my $idex=${$orderedexons->{$chr}{"id"}}[$i];
		my $strandex=${$orderedexons->{$chr}{"strand"}}[$i];
		
		my $j=$firstreg;
		
		while($j<$nbreg && ${$repeats->{$chr}{"end"}}[$j] < $startex){
		    $j++;
		}
		
		$firstreg=$j;
		
		while($j<$nbreg && ${$repeats->{$chr}{"start"}}[$j] <= $endex){
		    my $startreg=${$repeats->{$chr}{"start"}}[$j];
		    my $endreg=${$repeats->{$chr}{"end"}}[$j];
		    my $strandreg=${$repeats->{$chr}{"strand"}}[$j];
		    
		    my $M=max($startex, $startreg);
		    my $m=min($endex, $endreg);
		    
		    if($M<=$m){
			if(($type eq "sense" && $strandex eq $strandreg) || ($type eq "antisense" && $strandex ne $strandreg) || $type eq "any"){
			    if(exists $overlap->{$idex}){
				push(@{$overlap->{$idex}{"start"}}, $M);
				push(@{$overlap->{$idex}{"end"}}, $m);
			    }
			    else{
				$overlap->{$idex}={"start"=>[$M], "end"=>[$m]};
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
    print "This script computes overlap with repeats.\n";
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
$parameters{"pathRepeatMasker"}="NA";
$parameters{"pathOutputSense"}="NA";
$parameters{"pathOutputAntisense"}="NA";
$parameters{"pathOutputBothStrands"}="NA";

my @defaultpars=("pathAnnotGTF", "pathRepeatMasker", "pathOutputSense", "pathOutputAntisense", "pathOutputBothStrands");

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

my %repeats;

readRepeatMasker($parameters{"pathRepeatMasker"}, \%repeats);

print "Done.\n";

##############################################################

print "Ordering exons...\n";

my %orderedexons;

orderExons(\%exons, \%orderedexons);

my %orderedrepeats;

orderExons(\%repeats, \%orderedrepeats);

print "Done.\n";

##############################################################

print "Computing overlap with repeats...\n";

my %overlapsense;

computeOverlap(\%orderedexons, \%orderedrepeats, "sense", \%overlapsense);

my $nbov=keys %overlapsense;

print "Found ".$nbov." exons with repeat coverage on the sense strand.\n";


my %overlapantisense;

computeOverlap(\%orderedexons, \%orderedrepeats, "antisense", \%overlapantisense);

my $nbov=keys %overlapantisense;

print "Found ".$nbov." exons with repeat coverage on the antisense strand.\n";


my %overlapboth;

computeOverlap(\%orderedexons, \%orderedrepeats, "any", \%overlapboth);

my $nbov=keys %overlapboth;

print "Found ".$nbov." exons with repeat coverage on any strand.\n";


   
print "Done.\n";

##############################################################

print "Writing output...\n";

## sense

open(my $output, ">".$parameters{"pathOutputSense"});

print $output "GeneID\tTranscriptID\tExonID\tChr\tStart\tEnd\tStrand\n";

foreach my $txid (keys %txexons){
    foreach my $exonid (keys %{$txexons{$txid}}){
	my %blocks;
	
	my $geneid=join(";", keys %{$exons{$exonid}{"genes"}});
	my $chr=$exons{$exonid}{"chr"};
	my $start=$exons{$exonid}{"start"};
	my $end=$exons{$exonid}{"end"};
	my $strand=$exons{$exonid}{"strand"};
	
	my $totlen=$end-$start+1;
	
	if(exists $overlapsense{$exonid}){
	    makeBlocks($overlapsense{$exonid}, \%blocks);

	    my $nbex=@{$blocks{"start"}};
	    
	    for(my $i=0; $i<$nbex; $i++){
		print $output $geneid."\t".$txid."\t".$exonid."\t".$chr."\t".${$blocks{"start"}}[$i]."\t".${$blocks{"end"}}[$i]."\t".$strand."\n";
	    }
	}
    }
}

close($output);

## antisense

open(my $output, ">".$parameters{"pathOutputAntisense"});

print $output "GeneID\tTranscriptID\tExonID\tChr\tStart\tEnd\tStrand\n";

foreach my $txid (keys %txexons){
    foreach my $exonid (keys %{$txexons{$txid}}){
	my %blocks;
	
	my $geneid=join(";", keys %{$exons{$exonid}{"genes"}});
	my $chr=$exons{$exonid}{"chr"};
	my $start=$exons{$exonid}{"start"};
	my $end=$exons{$exonid}{"end"};
	my $strand=$exons{$exonid}{"strand"};
	
	my $totlen=$end-$start+1;
	
	if(exists $overlapantisense{$exonid}){
	    makeBlocks($overlapantisense{$exonid}, \%blocks);

	    my $nbex=@{$blocks{"start"}};
	    
	    for(my $i=0; $i<$nbex; $i++){
		print $output $geneid."\t".$txid."\t".$exonid."\t".$chr."\t".${$blocks{"start"}}[$i]."\t".${$blocks{"end"}}[$i]."\t".$strand."\n";
	    }
	}
    }
}

close($output);

## any strand

open(my $output, ">".$parameters{"pathOutputBothStrands"});

print $output "GeneID\tTranscriptID\tExonID\tChr\tStart\tEnd\tStrand\n";

foreach my $txid (keys %txexons){
    foreach my $exonid (keys %{$txexons{$txid}}){
	my %blocks;
	
	my $geneid=join(";", keys %{$exons{$exonid}{"genes"}});
	my $chr=$exons{$exonid}{"chr"};
	my $start=$exons{$exonid}{"start"};
	my $end=$exons{$exonid}{"end"};
	my $strand=$exons{$exonid}{"strand"};
	
	my $totlen=$end-$start+1;
	
	if(exists $overlapboth{$exonid}){
	    makeBlocks($overlapboth{$exonid}, \%blocks);

	    my $nbex=@{$blocks{"start"}};
	    
	    for(my $i=0; $i<$nbex; $i++){
		print $output $geneid."\t".$txid."\t".$exonid."\t".$chr."\t".${$blocks{"start"}}[$i]."\t".${$blocks{"end"}}[$i]."\t".$strand."\n";
	    }
	}
    }
}

close($output);

print "Done.\n";

##############################################################
