use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $exons=$_[1];
    my $txgene=$_[2];
    
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

	    my $info=$s[8];
	    my @infoarray=split(";", $info);

	    my $geneid=findInfo("gene_id", \@infoarray);
	    my $txid=findInfo("transcript_id", \@infoarray);

	    $txgene->{$txid}=$geneid;
	    
	    my $exonid=$chr.",".$start.",".$end.",".$strand;

	    if(exists $exons->{$exonid}){
		$exons->{$exonid}{"transcripts"}{$txid}=1;
		$exons->{$exonid}{"genes"}{$geneid}=1;
	    } else{
		$exons->{$exonid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand, "genes"=>{$geneid=>1}, "transcripts"=>{$txid=>1}};
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
    
    my @grepres=grep(/${pattern}/,@{$array});

    my $nbg=@grepres;
    my $nbreal=0;

    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    } else{
	my $nbreal=0;
	
	foreach my $g (@grepres){
	    $g =~ s/^\s+//; ## remove whitespace
	    my @u=split(" ",$g);

	    if($u[0] eq $pattern){
		$nbreal++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }
    
    if($nbreal>1){
	return "NA";
    }
    return $res;
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

sub readCoverage{
    my $pathin=$_[0];
    my $coverage=$_[1];

    my $input;
    my @s=split("\\.",$pathin);
    my $ext=$s[-1];
    
    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    }
    else{
	open($input, $pathin);
    }
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t" , $line);
	
	my $chr=$s[0];
	my $start=$s[1]+1; ## 0-based, included
	my $end=$s[2]; ## 0-based, not included
	my $score=$s[3]+0.0;
	
	if(exists $coverage->{$chr}){
	    my $laststart=${$coverage->{$chr}{"start"}}[-1];
	    
	    if($start<$laststart){
		print "Data are not ordered! ".$start." ".$laststart."\n";
		exit(1);
	    }
	    
	    push(@{$coverage->{$chr}{"start"}}, $start);
	    push(@{$coverage->{$chr}{"end"}}, $end);
	    push(@{$coverage->{$chr}{"score"}}, $score);
	}
	else{
	    $coverage->{$chr}={"start"=>[$start], "end"=>[$end], "score"=>[$score]};
	}
	
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub computeExonCoverage{
    my $orderedexons=$_[0];
    my $coverage=$_[1];
    my $covexons=$_[2];

    foreach my $chr (keys %{$orderedexons}){
	foreach my $strand (keys %{$orderedexons->{$chr}}){
	    print $chr."\t".$strand."\n";

	    if(exists $coverage->{$chr}){
		my $nbex=@{$orderedexons->{$chr}{$strand}{"start"}};
		my $nbreg=@{$coverage->{$chr}{"start"}};
		
		my $firstreg=0;
		
		for(my $i=0; $i<$nbex; $i++){
		    my $startex=${$orderedexons->{$chr}{$strand}{"start"}}[$i];
		    my $endex=${$orderedexons->{$chr}{$strand}{"end"}}[$i];
		    my $idex=${$orderedexons->{$chr}{$strand}{"id"}}[$i];
		    
		    my $lenex=($endex-$startex+1.0);
		    my $sumcov=0;
		    
		    my $j=$firstreg;
		    
		    while($j<$nbreg && ${$coverage->{$chr}{"end"}}[$j] < $startex){
			$j++;
		    }

		    $firstreg=$j;
		    
		    while($j<$nbreg && ${$coverage->{$chr}{"start"}}[$j] <= $endex){
			my $startreg=${$coverage->{$chr}{"start"}}[$j];
			my $endreg=${$coverage->{$chr}{"end"}}[$j];
			
			my $M=max($startex, $startreg);
			my $m=min($endex, $endreg);
			
			my $lenov=($m-$M+1);
			
			if($lenov>=1){
			    my $score=${$coverage->{$chr}{"score"}}[$j];
			    
			    $sumcov+=($lenov+0.0)*($score+0.0);
			}
			
			$j++;
		    }

		    my $meancov=($sumcov+0.0)/($lenex);
		    
		    $covexons->{$idex}=$meancov;
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
    print "This script computes exon RNA-seq coverage.\n";
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
$parameters{"pathCoverageForward"}="NA";
$parameters{"pathCoverageReverse"}="NA";
$parameters{"pathOutputExons"}="NA";
$parameters{"pathOutputTranscripts"}="NA";

my @defaultpars=("pathAnnotGTF", "pathCoverageForward", "pathCoverageReverse", "pathOutputExons", "pathOutputTranscripts");

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
my %txgene;
readGTF($parameters{"pathAnnotGTF"},  \%exons, \%txgene);

my $nbex=keys %exons;
my $nbtx=keys %txgene;

print "Found ".$nbex." exons and ".$nbtx." transcripts.\n";

print "Done.\n";


print "Ordering exons...\n";

my %orderedexons;

orderExons(\%exons, \%orderedexons);

print "Done.\n";

##############################################################

print "Reading RNA-seq coverage...\n";

my %coverageforward;
my %coveragereverse;

readCoverage($parameters{"pathCoverageForward"}, \%coverageforward);
readCoverage($parameters{"pathCoverageReverse"}, \%coveragereverse);

print "Done.\n";

##############################################################

print "Computing exon coverage...\n";

my %exoncoveragefwd;
my %exoncoveragerev;

print "forward\n";
computeExonCoverage(\%orderedexons, \%coverageforward, \%exoncoveragefwd);

print "reverse\n";
computeExonCoverage(\%orderedexons, \%coveragereverse, \%exoncoveragerev);

print "Done.\n";

##############################################################

print "Writing output for exons...\n";

open(my $output, ">".$parameters{"pathOutputExons"});

print $output "ExonID\tChr\tStart\tEnd\tStrand\tGeneID\tTranscriptID\tTotalLength\tCoverageSense\tCoverageAntisense\n";

my %txsense;
my %txantisense;

foreach my $exonid (keys %exons){
    my $geneid=join(";", keys %{$exons{$exonid}{"genes"}});
    my $txid=join(";", keys %{$exons{$exonid}{"transcripts"}});
    my $chr=$exons{$exonid}{"chr"};
    my $start=$exons{$exonid}{"start"};
    my $end=$exons{$exonid}{"end"};
    my $strand=$exons{$exonid}{"strand"};

    my $totlen=$end-$start+1;
    my $covsense=0;
    my $covantisense=0;
    
    if($strand eq "1"){
	if(exists $exoncoveragefwd{$exonid}){
	    $covsense=$exoncoveragefwd{$exonid};
	}

	if(exists $exoncoveragerev{$exonid}){
	    $covantisense=$exoncoveragerev{$exonid};
	}
    }

    if($strand eq "-1"){
	if(exists $exoncoveragefwd{$exonid}){
	    $covantisense=$exoncoveragefwd{$exonid};
	}

	if(exists $exoncoveragerev{$exonid}){
	    $covsense=$exoncoveragerev{$exonid};
	}
    }

    print $output $exonid."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$geneid."\t".$txid."\t".$totlen."\t".$covsense."\t".$covantisense."\n";

    foreach my $txid (keys %{$exons{$exonid}{"transcripts"}}){
	if(exists $txsense{$txid}){
	    $txsense{$txid}{"sumcov"}+=($covsense+0.0)*($totlen+0.0);
	    $txsense{$txid}{"length"}+=$totlen;
	}
	else{
	    $txsense{$txid}={"sumcov"=>($covsense+0.0)*($totlen+0.0), "length"=>$totlen};
	}

	if(exists $txantisense{$txid}){
	    $txantisense{$txid}{"sumcov"}+=($covantisense+0.0)*($totlen+0.0);
	    $txantisense{$txid}{"length"}+=$totlen;
	}
	else{
	    $txantisense{$txid}={"sumcov"=>($covantisense+0.0)*($totlen+0.0), "length"=>$totlen};
	}
    }
}

close($output);


my $nbs=keys %txsense;
my $nba=keys %txantisense;

print "Found ".$nbs." transcripts with sense coverage and ".$nba." transcripts with antisense coverage.\n";
print "Done.\n";

##############################################################

print "Writing output for transcripts...\n";

open(my $output, ">".$parameters{"pathOutputTranscripts"});

print $output "TranscriptID\tGeneID\tTotalLength\tCoverageSense\tCoverageAntisense\n";

foreach my $txid (keys %txsense){
    if(!exists $txgene{$txid}){
	print "Weird! cannot find gene info for ".$txid."\n";
	exit(1);
    }
    
    my $geneid=$txgene{$txid};

    my $totlen1=$txsense{$txid}{"length"};
    my $totlen2=$txantisense{$txid}{"length"};

    if($totlen1!=$totlen2){
	print "Weird! different sense/antisense length for ".$txid."\n";
	exit(1);
    }
    
    my $covsense=($txsense{$txid}{"sumcov"}+0.0)/($txsense{$txid}{"length"}+0.0);
    my $covantisense=($txantisense{$txid}{"sumcov"}+0.0)/($txantisense{$txid}{"length"}+0.0);

    print $output $txid."\t".$geneid."\t".$totlen1."\t".$covsense."\t".$covantisense."\n";
}

print "Done.\n";

##############################################################
