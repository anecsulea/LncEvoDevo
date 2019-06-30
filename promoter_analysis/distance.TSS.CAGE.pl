#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $exoncoords=$_[1];
    my $genetx=$_[2];
    my $txex=$_[3];
  
    open(my $input, $pathin);
    
    my $line=<$input>;
    my $prefix=substr $line, 0,1;
    
    while($prefix eq "#"){
	$line=<$input>;
	$prefix=substr $line, 0,1;
    }

    my $nbunstranded=0; ## we discard genes with "." strand

    my %rttx;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon"){
	
	    my $chr=$s[0];
	    my $start=$s[3]+0; ## 1-based
	    my $end=$s[4]+0;
	    my $strand=$s[6];
	  
	    if($strand ne "."){
		my $info=$s[8];
		my @t=split(";", $info);
		my $geneid=findInfo("gene_id", \@t);
		my $txid=findInfo("transcript_id", \@t);
		
		if($txid eq "NA"){
		    print "could not find transcript in ".$line."\n";
		    exit(1);
		}		

		if($geneid eq "NA"){
		    print "could not find gene in ".$line."\n";
		    exit(1);
		}		
		
		my $exonid=$chr.",".$start.",".$end.",".$strand;
		
		## fill in exon coords
		
		$exoncoords->{$exonid}={"chr"=>$chr,"start"=>$start, "end"=>$end, "strand"=>$strand};
		
		## exon - tx correspondence
	
		if(exists $genetx->{$geneid}){
		    $genetx->{$geneid}{$txid}=1;
		}
		else{
		    $genetx->{$geneid}={$txid=>1};
		}
		
		if(exists $txex->{$txid}){
		    $txex->{$txid}{$exonid}=1;
		}
		else{
		    $txex->{$txid}={$exonid=>1};
		}
	    }
	    else{
		$nbunstranded++;
	    }
	}
   	
	$line=<$input>;
    }
    
    close($input);

    print "We discarded ".$nbunstranded." exons with undefined strand.\n";
  
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

##########################################################################

sub extractTSS{
    my $txexons=$_[0];
    my $exoncoords=$_[1];
    my $tsscoords=$_[2];

    foreach my $tx (keys %{$txexons}){
	my $chr="NA";
	my $start="NA";
	my $end="NA";
	my $strand="NA";

	foreach my $exid (keys %{$txexons->{$tx}}){
	    my $thischr=$exoncoords->{$exid}{"chr"};
	    my $thisstart=$exoncoords->{$exid}{"start"};
	    my $thisend=$exoncoords->{$exid}{"end"};
	    my $thisstrand=$exoncoords->{$exid}{"strand"};

	    if($chr ne "NA" && $chr ne $thischr){
		print "Weird! different chromosomes for ".$tx."\n";
		exit(1);
	    }

	    if($strand ne "NA" && $strand ne $thisstrand){
		print "Weird! different strands for ".$tx."\n";
		exit(1);
	    }
	    
	    $chr=$thischr;
	    $strand=$thisstrand;

	    if($start ne "NA"){
		if($thisstart<$start){
		    $start=$thisstart;
		}
	    } else{
		$start=$thisstart;
	    }

	    if($end ne "NA"){
		if($thisend>$end){
		    $end=$thisend;
		}
	    } else{
		$end=$thisend;
	    }
	}
	
	$tsscoords->{$tx}={"chr"=>$chr, "strand"=>$strand};
	
	if($strand eq "+" || $strand eq "1"){
	    $tsscoords->{$tx}{"position"}=$start;
	} else{
	    if($strand eq "-" || $strand eq "-1"){
		$tsscoords->{$tx}{"position"}=$end;
	    } else{
		print "Weird strand for ".$tx."\n";
	    }
	}
    }
}

##########################################################################

sub readCAGE{
    my $pathin=$_[0];
    my $cage=$_[1];

    my $input;
    my @s=split("\\.",$pathin);
    my $ext=$s[-1];

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
	    $chr=substr $chr, 3;
	}

	my $start=$s[1]+1;
	my $end=$s[2]+0;
	my $strand=$s[5];

	my $id=$chr.",".$start.",".$end.",".$strand;

	$cage->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};

	$line=<$input>;

    }
    
    close($input);
}

##########################################################################

sub combineCoordinates{
    my $tsscoords=$_[0];
    my $cagecoords=$_[1];
    my $combined=$_[2];

    my %unordered;

    foreach my $tx (keys %{$tsscoords}){
	my $chr=$tsscoords->{$tx}{"chr"};
	my $pos=$tsscoords->{$tx}{"position"};
	my $strand=$tsscoords->{$tx}{"strand"};

	if(exists $unordered{$chr}){
	    if(exists $unordered{$chr}{$strand}){
		if(exists $unordered{$chr}{$strand}{$pos}){
		    push(@{$unordered{$chr}{$strand}{$pos}{"id"}}, $tx);
		    push(@{$unordered{$chr}{$strand}{$pos}{"type"}}, "TSS");
		} else{
		    $unordered{$chr}{$strand}{$pos}={"id"=>[$tx], "type"=>["TSS"]};
		}
	    } else{
		$unordered{$chr}{$strand}={$pos=>{"id"=>[$tx], "type"=>["TSS"]}};
	    }
	} else{
	    $unordered{$chr}={$strand=>{$pos=>{"id"=>[$tx], "type"=>["TSS"]}}};
	}
    }
    
    foreach my $id (keys %{$cagecoords}){
	my $chr=$cagecoords->{$id}{"chr"};
	my $start=$cagecoords->{$id}{"start"};
	my $end=$cagecoords->{$id}{"end"};
	my $strand=$cagecoords->{$id}{"strand"};

	if(exists $unordered{$chr}){
	    if(exists $unordered{$chr}{$strand}){
		if(exists $unordered{$chr}{$strand}{$start}){
		    push(@{$unordered{$chr}{$strand}{$start}{"id"}}, $id);
		    push(@{$unordered{$chr}{$strand}{$start}{"type"}}, "CAGE");
		} else{
		    $unordered{$chr}{$strand}{$start}={"id"=>[$id], "type"=>["CAGE"]};
		}

		if(exists $unordered{$chr}{$strand}{$end}){
		    push(@{$unordered{$chr}{$strand}{$end}{"id"}}, $id);
		    push(@{$unordered{$chr}{$strand}{$end}{"type"}}, "CAGE");
		} else{
		    $unordered{$chr}{$strand}{$end}={"id"=>[$id], "type"=>["CAGE"]};
		}
		
	    } else{
		$unordered{$chr}{$strand}={$start=>{"id"=>[$id], "type"=>["CAGE"]}, $end=>{"id"=>[$id], "type"=>["CAGE"]}};
	    }
	} else{
	    $unordered{$chr}={$strand=>{$start=>{"id"=>[$id], "type"=>["CAGE"]}, $end=>{"id"=>[$id], "type"=>["CAGE"]}}};
	}
    }

    ## now order positions


    foreach my $chr (keys %unordered){
	$combined->{$chr}={};
	
	foreach my $strand (keys %{$unordered{$chr}}){
	    $combined->{$chr}{$strand}={"id"=>[], "position"=>[], "type"=>[]};

	    my @unsortedpos=keys %{$unordered{$chr}{$strand}};
	    my @sortedpos=sort {$a<=>$b} @unsortedpos;

	    foreach my $pos (@sortedpos){
		my $nb=@{$unordered{$chr}{$strand}{$pos}{"id"}};

		for(my $i=0; $i<$nb; $i++){
		    push(@{$combined->{$chr}{$strand}{"position"}}, $pos);
		    push(@{$combined->{$chr}{$strand}{"id"}}, ${$unordered{$chr}{$strand}{$pos}{"id"}}[$i]);
		    push(@{$combined->{$chr}{$strand}{"type"}}, ${$unordered{$chr}{$strand}{$pos}{"type"}}[$i]);
		}
	    }
	}
    }
}

##########################################################################

sub computeDistanceTSSCAGE{
    my $combined=$_[0];
    my $distance=$_[1];

    foreach my $chr (keys %{$combined}){
	foreach my $strand (keys %{$combined->{$chr}}){
	    my $nb=@{$combined->{$chr}{$strand}{"position"}};

	    for(my $i=0; $i<$nb; $i++){
		my $thispos=${$combined->{$chr}{$strand}{"position"}}[$i];
		my $thisid=${$combined->{$chr}{$strand}{"id"}}[$i];
		my $thistype=${$combined->{$chr}{$strand}{"type"}}[$i];

		if($thistype eq "TSS"){
		    my $distleft="Inf";
		    my $distright="Inf";
		    my $idleft="Inf";
		    my $idright="Inf";

		    for(my $j=($i-1); $j>=0; $j--){
			my $typeleft=${$combined->{$chr}{$strand}{"type"}}[$j];

			if($typeleft eq "CAGE"){
			    $idleft=${$combined->{$chr}{$strand}{"id"}}[$j];
			    my $posleft=${$combined->{$chr}{$strand}{"position"}}[$j];
			    $distleft=$thispos-$posleft;
			    
			    last;
			}
		    }

		    for(my $j=($i+1); $j<$nb; $j++){
			my $typeright=${$combined->{$chr}{$strand}{"type"}}[$j];

			if($typeright eq "CAGE"){
			    $idright=${$combined->{$chr}{$strand}{"id"}}[$j];
			    my $posright=${$combined->{$chr}{$strand}{"position"}}[$j];
			    $distright=$posright-$thispos;
			    
			    last;
			}
		    }

		    if($distleft ne "Inf"){
			if($distright ne "Inf"){
			    if($distleft<=$distright){
				$distance->{$thisid}={"distance"=>$distleft,"closest"=>$idleft};
			    } else{
				$distance->{$thisid}={"distance"=>$distright,"closest"=>$idright};
			    }
			} else{
			    $distance->{$thisid}={"distance"=>$distleft,"closest"=>$idleft};
			}
		    } else{
			if($distright ne "Inf"){
			    $distance->{$thisid}={"distance"=>$distright,"closest"=>$idright};
			} else{
			    $distance->{$thisid}={"distance"=>"Inf","closest"=>"NA"};
			}
		    }
		}
	    }
	}
    }
}

##########################################################################

sub orderTSS{
    my $unordered=$_[0];
    my $ordered=$_[1];

    my %hashpos;
    
    foreach my $tx (keys %{$unordered}){
	my $chr=$unordered->{$tx}{"chr"};
	my $strand=$unordered->{$tx}{"strand"};
	my $pos=$unordered->{$tx}{"position"};
	
	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$strand}){
		if(exists $hashpos{$chr}{$strand}{$pos}){
		    push(@{$hashpos{$chr}{$strand}{$pos}}, $tx);
		} else{
		    $hashpos{$chr}{$strand}{$pos}=[$tx];
		}
	    } else{
		$hashpos{$chr}{$strand}={$pos=>[$tx]};
	    }
	} else{
	    $hashpos{$chr}={$strand=>{$pos=>[$tx]}};
	}
    }

    foreach my $chr (keys %hashpos){
	$ordered->{$chr}={};

	foreach my $strand (keys %{$hashpos{$chr}}){
	    $ordered->{$chr}{$strand}={"position"=>[], "id"=>[]};

	    my @unsortedpos=keys %{$hashpos{$chr}{$strand}};
	    my @sortedpos=sort {$a<=>$b} @unsortedpos;

	    foreach my $pos (@sortedpos){
		my $nb=@{$hashpos{$chr}{$strand}{$pos}};

		for(my $i=0; $i<$nb; $i++){
		    push(@{$ordered->{$chr}{$strand}{"position"}}, $pos);
		    push(@{$ordered->{$chr}{$strand}{"id"}}, ${$hashpos{$chr}{$strand}{$pos}}[$i]);
		}
	    }
	}
    }
}

##########################################################################

sub orderCAGE{
    my $unordered=$_[0];
    my $ordered=$_[1];

    my %hashpos;
    
    foreach my $id (keys %{$unordered}){
	my $chr=$unordered->{$id}{"chr"};
	my $strand=$unordered->{$id}{"strand"};
	my $start=$unordered->{$id}{"start"};
	my $end=$unordered->{$id}{"end"};
	
	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$strand}){
		if(exists $hashpos{$chr}{$strand}{$start}){
		    push(@{$hashpos{$chr}{$strand}{$start}["id"]}, $id);
		    push(@{$hashpos{$chr}{$strand}{$start}["end"]}, $end);
		} else{
		    $hashpos{$chr}{$strand}{$start}={"id"=>[$id], "end"=>[$end]};
		}
	    } else{
		$hashpos{$chr}{$strand}={$start=>{"id"=>[$id], "end"=>[$end]}};
	    }
	} else{
	    $hashpos{$chr}={$strand=>{$start=>{"id"=>[$id], "end"=>[$end]}}};
	}
    }

    foreach my $chr (keys %hashpos){
	$ordered->{$chr}={};

	foreach my $strand (keys %{$hashpos{$chr}}){
	    $ordered->{$chr}{$strand}={"position"=>[], "id"=>[]};

	    my @unsortedpos=keys %{$hashpos{$chr}{$strand}};
	    my @sortedpos=sort {$a<=>$b} @unsortedpos;

	    foreach my $pos (@sortedpos){
		my $nb=@{$hashpos{$chr}{$strand}{$pos}{"id"}};

		for(my $i=0; $i<$nb; $i++){
		    push(@{$ordered->{$chr}{$strand}{"start"}}, $pos);
		    push(@{$ordered->{$chr}{$strand}{"id"}}, ${$hashpos{$chr}{$strand}{$pos}{"id"}}[$i]);
		    push(@{$ordered->{$chr}{$strand}{"end"}}, ${$hashpos{$chr}{$strand}{$pos}{"end"}}[$i]);
		}
	    }
	}
    }
}

##########################################################################

sub extractOverlap{
    my $reftss=$_[0]; ## ordered TSS
    my $refcage=$_[1]; ## ordered CAGE peaks
    my $refoverlap=$_[2]; ## overlap with blocks
   
    foreach my $chr (keys %{$reftss}){
	foreach my $strand (keys %{$reftss->{$chr}}){
	    my $nbtss=@{$reftss->{$chr}{$strand}{"position"}};
	    
	    if(exists $refcage->{$chr} && exists $refcage->{$chr}{$strand}){
		my $nbcage=@{$refcage->{$chr}{$strand}{"start"}};
		
		my $firstindex=0;  ## this is where we start looking for overlap
	    
		for(my $i=0; $i<$nbtss; $i++){
		    my $pos=${$reftss->{$chr}{$strand}{"position"}}[$i];
		    my $idtx=${$reftss->{$chr}{$strand}{"id"}}[$i];
		   		    
		    my $j=$firstindex;
		    
		    while($j<$nbcage && ${$refcage->{$chr}{$strand}{"end"}}[$j]<$pos){ ## there cannnot be any overlap before that 
			$j++;
		    }
		    
		    $firstindex=$j;
		    
		    while($j<$nbcage && ${$refcage->{$chr}{$strand}{"start"}}[$j]<=$pos){  ## we stop looking for overlap if the start coordinate of the second set of blocks is larger than the end coordinate of block
			my $idcage=${$refcage->{$chr}{$strand}{"id"}}[$j];
			
			if(exists $refoverlap->{$idtx}){
			    push(@{$refoverlap->{$idtx}}, $idcage);
			}
			else{
			    $refoverlap->{$idtx}=[$idcage];
			}
			
			$j++;
		    }
		}
	    }
	}
    }
}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes smallest distance between TSS and CAGE peaks.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################
##########################################################################

## parameters 

my %parameters;
$parameters{"pathGTF"}="NA";
$parameters{"pathCAGE"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathGTF", "pathCAGE", "pathOutput");


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

my %exoncoords;
my %genetx;
my %txexons;

readGTF($parameters{"pathGTF"}, \%exoncoords, \%genetx, \%txexons);
my $nbg=keys %genetx;
my $nbex=keys %exoncoords;
my $nbtx=keys %txexons;

print "Found ".$nbg." genes, ".$nbtx." transcripts and ".$nbex." exons.\n";
    
print "Done.\n";

##############################################################

print "Extracting and ordering TSS coordinates...\n";

my %tsscoords;
extractTSS(\%txexons, \%exoncoords, \%tsscoords);
my $nbtss=keys %tsscoords;

print "Found ".$nbtss." TSS.\n";

my %orderedtss;
orderTSS(\%tsscoords, \%orderedtss);

print "Done.\n";

##############################################################

print "Reading and ordering CAGE peaks...\n";

my %cage;
readCAGE($parameters{"pathCAGE"}, \%cage);

my $nbc=keys %cage;
print "Found ".$nbc." CAGE peaks.\n";

my %orderedcage;
orderCAGE(\%cage, \%orderedcage);

print "Done.\n";

##############################################################

print "Extracting overlap between TSS and CAGE...\n";

my %overlaptsscage;

extractOverlap(\%orderedtss, \%orderedcage, \%overlaptsscage); 

print "Done.\n";

##############################################################

print "Combining and ordering coordinates...\n";

my %combined;
combineCoordinates(\%tsscoords, \%cage, \%combined);
 
print "Done.\n";

##############################################################

print "Computing distance between TSS and CAGE peaks...\n";

my %distcage;
computeDistanceTSSCAGE(\%combined, \%distcage);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});
print $output "GeneID\tTranscriptID\tChr\tStrand\tTSS\tDistance\tIDClosestPeak\n";

foreach my $gene (keys %genetx){
    foreach my $tx (keys %{$genetx{$gene}}){

	if(!exists $distcage{$tx}){
	    print "Weird! undefined distance for ".$tx."\n";
	    exit(1);
	}

	if(exists $overlaptsscage{$tx}){
	    print $output $gene."\t".$tx."\t".$tsscoords{$tx}{"chr"}."\t".$tsscoords{$tx}{"strand"}."\t".$tsscoords{$tx}{"position"}."\t-1\t".join(";",@{$overlaptsscage{$tx}})."\n";
	} else{
	
	    print $output $gene."\t".$tx."\t".$tsscoords{$tx}{"chr"}."\t".$tsscoords{$tx}{"strand"}."\t".$tsscoords{$tx}{"position"}."\t".$distcage{$tx}{"distance"}."\t".$distcage{$tx}{"closest"}."\n";
	}
    }
}


close($output);

print "Done.\n";

##############################################################
